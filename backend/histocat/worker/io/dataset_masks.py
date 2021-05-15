import logging
import os
from pathlib import Path
from typing import Dict, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from skimage import measure
from sqlalchemy.orm import Session
from tifffile import tifffile

from histocat.core.acquisition import service as acquisition_service
from histocat.core.acquisition.models import AcquisitionModel
from histocat.core.dataset import service as dataset_service
from histocat.core.dataset.dto import DatasetCreateDto, DatasetUpdateDto
from histocat.core.dataset.models import CELL_FILENAME, DatasetModel
from histocat.core.notifier import Message
from histocat.core.project import service as project_service
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.worker.io.utils import copy_file

logger = logging.getLogger(__name__)

MASKS_CSV_FILE = "masks.csv"


def import_dataset(db: Session, mask_map_filename: Path, project_id: int):
    """Import dataset from the segmentation masks only."""

    project = project_service.get(db, id=project_id)
    if not project:
        logger.warning(f"Cannot import dataset: project [id: {project_id}] does not exist.")
        return

    create_params = DatasetCreateDto(
        project_id=project_id, name=mask_map_filename.stem, origin="Mask", status="pending"
    )
    dataset = dataset_service.create(db, params=create_params)

    src_folder = mask_map_filename.parent
    dst_folder = Path(dataset.location)
    os.makedirs(dst_folder, exist_ok=True)

    mask_map = _import_mask_map_csv(mask_map_filename)
    masks = {}
    segmentation_data = []
    for index, row in mask_map.iterrows():
        acquisition, mask_meta = _import_mask(db, project_id, src_folder, dst_folder, row["Image"], row["Mask"])
        if mask_meta is not None:
            masks[acquisition.id] = mask_meta
            ac_output = _process_single_cell_data(acquisition, mask_meta.get("location"))
            segmentation_data.append(ac_output)

    meta = {"masks": masks, "columns": {}}

    channels = _import_cells(dataset=dataset, segmentation_data=segmentation_data)

    acquisition_ids = sorted(list(masks.keys()))

    update_params = DatasetUpdateDto(status="ready", acquisition_ids=acquisition_ids, channels=channels, meta=meta)
    dataset = dataset_service.update(db, item=dataset, params=update_params)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(project_id, "dataset_imported"))


def _import_mask_map_csv(mask_map_filename: Path):
    df = pd.read_csv(mask_map_filename)
    return df


def _import_mask(
    db: Session, project_id: int, src_folder: Path, dst_folder: Path, acquisition_name: str, mask_filename: str
):
    acquisition = acquisition_service.get_by_project_id_and_description(
        db, project_id=project_id, description=acquisition_name
    )
    uri = src_folder / mask_filename
    location = copy_file(str(uri), str(dst_folder))
    mask_meta = {
        "location": location,
        "slide": {"id": acquisition.slide.id, "origin_id": acquisition.slide.origin_id},
        "acquisition": {"id": acquisition.id, "origin_id": acquisition.origin_id},
    }
    return acquisition, mask_meta


def _process_single_cell_data(acquisition: AcquisitionModel, mask_location: str):
    ometiff = OmeTiffParser(acquisition.location)
    mask = tifffile.imread(mask_location)
    acquisition_data = ometiff.get_acquisition_data()

    output = {
        "acquisition": acquisition,
        "mask_location": mask_location,
        "object_numbers": None,
        "centroids_x": None,
        "centroids_y": None,
        "mean_intensities": {},
        "channel_names": acquisition_data.channel_names,
    }

    for c in acquisition_data.channel_names:
        d = measure.regionprops_table(
            label_image=mask,
            intensity_image=acquisition_data.get_image_by_name(c),
            properties=("label", "centroid", "mean_intensity"),
        )
        if output["object_numbers"] is None:
            output["object_numbers"] = d.get("label")
            output["centroids_x"] = d.get("centroid-1")
            output["centroids_y"] = d.get("centroid-0")
        output["mean_intensities"][c] = d.get("mean_intensity")

    return output


def _import_cells(dataset: DatasetModel, segmentation_data: Sequence[Dict]):
    dst = os.path.join(dataset.location, CELL_FILENAME)

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
    obs.index = obs.index.astype(str, copy=False)
    obs["CellId"] = obs.index

    var_names = []
    x_df = pd.DataFrame()
    for k, v in mean_intensities_all.items():
        x_df[k] = v
        var_names.append(k)
    var = pd.DataFrame(index=var_names)
    var.index.astype(str, copy=False)
    var["Channel"] = var.index
    X_counts = x_df.to_numpy()

    adata = ad.AnnData(X_counts, obs=obs, var=var, dtype="float32")
    adata.write_h5ad(dst)

    return var_names
