import logging
import os
from typing import Sequence

import cv2
import numpy as np
import scanpy as sc
from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import ORJSONResponse
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from skimage.measure import regionprops
from sqlalchemy.orm import Session
from starlette.requests import Request
from sklearn.ensemble import RandomForestClassifier
from starlette import status

from histocat.api.db import get_db
from histocat.api.security import get_active_member
from histocat.core.acquisition import service as acquisition_service
from histocat.core.constants import ANNDATA_FILE_EXTENSION
from histocat.core.dataset import service as dataset_service
from histocat.core.result import service as result_service
from histocat.core.analysis.dto import RegionChannelStatsDto, RegionStatsSubmissionDto, ClassifyCellsSubmissionDto, \
    ClassifyCellsDto
from histocat.core.member.models import MemberModel

logger = logging.getLogger(__name__)

router = APIRouter()


@router.post("/groups/{group_id}/analysis/region", response_model=Sequence[RegionChannelStatsDto])
async def calculate_region_stats(
    params: RegionStatsSubmissionDto,
    request: Request,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """
    Calculate region's statistics
    """

    acquisition = acquisition_service.get_by_id(db, params.acquisition_id)
    parser = OmeTiffParser(acquisition.location)
    acq = parser.get_acquisition_data()
    mask = None
    contour = np.array(params.region_polygon).astype(int)
    content = []
    for metal in acq.channel_names:
        channel_img = acq.get_image_by_name(metal)
        if mask is None:
            mask = np.zeros(channel_img.shape, np.uint8)
            mask = cv2.drawContours(mask, [contour], 0, 255, -1)
        props = regionprops(mask, intensity_image=channel_img, cache=True, coordinates=None)
        props = props[0]
        content.append(
            {
                "metal": metal,
                "min": float("{0:.3f}".format(props.min_intensity)),
                "max": float("{0:.3f}".format(props.max_intensity)),
                "mean": float("{0:.3f}".format(props.mean_intensity)),
            }
        )

    return ORJSONResponse(content)


@router.post("/groups/{group_id}/analysis/classify", response_model=ClassifyCellsDto)
async def classify_cells(
    group_id: int,
    params: ClassifyCellsSubmissionDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """Classify cells."""

    if params.result_id:
        result = result_service.get(db, id=params.result_id)
        location = os.path.join(result.location, f"output{ANNDATA_FILE_EXTENSION}")
        adata = sc.read_h5ad(location)
    else:
        dataset = dataset_service.get(db, id=params.dataset_id)
        adata = sc.read_h5ad(dataset.cell_file_location())

    cell_ids = []
    cell_classes = []
    for annotation in params.annotations:
        ann_cell_ids = annotation["cellIds"]
        cell_ids.extend(ann_cell_ids)
        ann_cell_classes = [annotation["cellClass"]] * len(ann_cell_ids)
        cell_classes.extend(ann_cell_classes)

    # Convert cell ids to strings
    cell_ids = list(map(str, cell_ids))

    # Map cell ids to cell classes
    cells = dict(zip(cell_ids, cell_classes))

    df = adata.to_df()

    df_train = df[df.index.isin(cell_ids)].copy()
    if df_train.empty:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Train dataset is empty",
        )

    df_test = df[~df.index.isin(cell_ids)].copy()
    if df_test.empty:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Test dataset is empty",
        )

    df_train["cellClass"] = df_train.index
    df_train["cellClass"].replace(cells, inplace=True)

    # Create a Gaussian Classifier
    clf = RandomForestClassifier(n_estimators=params.n_estimators)

    # Train the model using the training sets y_pred=clf.predict(X_test)
    clf.fit(df_train[params.channels], df_train["cellClass"])

    y_pred = clf.predict(df_test[params.channels])
    y_pred_proba = clf.predict_proba(df_test[params.channels])

    df_test["cellClass"] = y_pred

    for index, cl in enumerate(clf.classes_):
        df_test[cl + "Prob"] = [prob[index] for prob in y_pred_proba]

    print(df_test)

    # Combine train and test dataframes together
    result_df = df_test.append(df_train)

    annotations = []
    for cell_class in result_df["cellClass"].unique():
        cellIds = result_df[result_df["cellClass"] == cell_class].index.to_list()
        annotation = {
            "cellClass": cell_class,
            "visible": True,
            "cellIds": cellIds
        }
        annotations.append(annotation)

    content = {
        "cellClasses": params.cell_classes,
        "annotations": annotations,
    }
    return ORJSONResponse(content)
