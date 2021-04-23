import logging
import os
from typing import Sequence

import cv2
import numpy as np
import scanpy as sc
from fastapi import APIRouter, Depends
from fastapi.responses import ORJSONResponse
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from skimage.measure import regionprops
from sqlalchemy.orm import Session
from starlette.requests import Request
from sklearn.ensemble import RandomForestClassifier

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

    cells = dict(zip(cell_ids, cell_classes))
    adata = adata[adata.obs["CellId"].isin(cell_ids)]
    df = adata.to_df()
    df["cellClass"] = df.index
    df["cellClass"] = df["cellClass"].astype(int)
    df["cellClass"].replace(cells, inplace=True)
    print(df)

    # Create a Gaussian Classifier
    # clf = RandomForestClassifier(n_estimators=100)

    # Train the model using the training sets y_pred=clf.predict(X_test)
    # clf.fit(X_train, y_train)
    #
    # y_pred = clf.predict(X_test)

    content = []
    return ORJSONResponse(content)
