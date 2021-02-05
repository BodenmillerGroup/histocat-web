from __future__ import annotations

import logging
import os
from typing import Dict, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from sqlalchemy.orm import Session

from histocat.core.acquisition.models import AcquisitionModel
from histocat.core.constants import ANNDATA_FILE_EXTENSION
from histocat.core.dataset.models import DatasetModel

logger = logging.getLogger(__name__)


EXPERIMENT_CSV_FILENAME = "Experiment"
ACQUISITION_METADATA_FILENAME = "acquisition_metadata"
OBJECT_RELATIONSHIPS_FILENAME = "Object relationships"
IMAGE_FILENAME = "Image"
CELL_FILENAME = "cell"

CHANNELS_FULL_CSV_ENDING = "_ac_full.csv"
PROBABILITIES_MASK_TIFF_ENDING = "_Probabilities_mask.tiff"


def import_dataset(db: Session, dataset: DatasetModel, segmentation_data: Sequence[Dict]):
    """Import dataset from the segmentation pipeline output."""

    meta = {"columns": {"neighbors": []}}

    masks = {}
    for ac_segmentation_data in segmentation_data:
        acquisition: AcquisitionModel = ac_segmentation_data.get("acquisition")
        mask_location = ac_segmentation_data.get("mask_location")
        masks[acquisition.id] = {
            "location": mask_location,
            "acquisition": {"id": acquisition.id, "origin_id": acquisition.origin_id},
            "slide": {"id": acquisition.slide.id, "origin_id": acquisition.slide.origin_id},
        }
    meta["masks"] = masks

    _import_cells(dataset=dataset, segmentation_data=segmentation_data)

    return meta


def _import_cells(dataset: DatasetModel, segmentation_data: Sequence[Dict]):
    dst = os.path.join(dataset.location, f"{CELL_FILENAME}{ANNDATA_FILE_EXTENSION}")

    object_numbers_all = None
    centroids_x_all = None
    centroids_y_all = None
    acquisition_ids_all = None
    image_numbers_all = None
    mean_intensities_all = {}

    for i, val in enumerate(segmentation_data):
        acquisition: AcquisitionModel = val.get("acquisition")

        object_numbers = val.get("object_numbers").astype(np.uint64)
        centroids_x = val.get("centroids_x").astype(np.float32)
        centroids_y = val.get("centroids_y").astype(np.float32)
        acquisition_ids = np.full(len(object_numbers), acquisition.id, np.uint64)
        image_numbers = np.full(len(object_numbers), i, np.uint64)

        if object_numbers_all is None:
            object_numbers_all = object_numbers
            centroids_x_all = centroids_x
            centroids_y_all = centroids_y
            acquisition_ids_all = acquisition_ids
            image_numbers_all = image_numbers
            for channel_name in val.get("channel_names"):
                mean_intensities_all[channel_name] = val.get("mean_intensities").get(channel_name).astype(np.float32)
        else:
            object_numbers_all = np.hstack((object_numbers_all, object_numbers))
            centroids_x_all = np.hstack((centroids_x_all, centroids_x))
            centroids_y_all = np.hstack((centroids_y_all, centroids_y))
            acquisition_ids_all = np.hstack((acquisition_ids_all, acquisition_ids))
            image_numbers_all = np.hstack((image_numbers_all, image_numbers))
            for channel_name in val.get("channel_names"):
                mean_intensities_all[channel_name] = np.hstack(
                    (
                        mean_intensities_all[channel_name],
                        val.get("mean_intensities").get(channel_name).astype(np.float32),
                    )
                )

    obs = pd.DataFrame(
        {
            "ObjectNumber": object_numbers_all,
            "AcquisitionId": acquisition_ids_all,
            "ImageNumber": image_numbers_all,
            "CentroidX": centroids_x_all,
            "CentroidY": centroids_y_all,
        }
    )
    obs["CellId"] = obs.index

    var_names = []
    x_df = pd.DataFrame()
    for k, v in mean_intensities_all.items():
        x_df[k] = v
        var_names.append(k)
    var = pd.DataFrame(index=var_names)
    var["Channel"] = var.index
    X_counts = x_df.to_numpy()

    adata = ad.AnnData(X_counts, obs=obs, var=var, dtype="float32")
    adata.write_h5ad(dst)

    return dst
