import os
from typing import Sequence, Union

import numpy as np
import tifffile
from deepcell.applications import Mesmer
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from skimage import measure
from sqlalchemy.orm import Session

from histocat.core.acquisition import service as acquisition_service
from histocat.core.dataset.models import DatasetModel
from histocat.core.errors import SegmentationError
from histocat.core.segmentation.dto import SegmentationSubmissionDto
from histocat.core.utils import timeit


def normalize_by_minmax(img: Sequence[Union[np.ndarray, np.ndarray]]):
    channel_mins = np.nanmin(img, axis=(1, 2), keepdims=True)
    channel_maxs = np.nanmax(img, axis=(1, 2), keepdims=True)
    img = (img - channel_mins) / (channel_maxs - channel_mins)
    return img


def normalize_by_zscore(img: Sequence[Union[np.ndarray, np.ndarray]]):
    channel_means = np.nanmean(img, axis=(1, 2), keepdims=True)
    channel_stds = np.nanstd(img, axis=(1, 2), keepdims=True)
    img = (img - channel_means) / channel_stds
    return img


@timeit
def process_acquisition(
    db: Session, acquisition_id: int, params: SegmentationSubmissionDto, model, dataset: DatasetModel
):
    acquisition = acquisition_service.get_by_id(db, acquisition_id)
    if not acquisition:
        raise SegmentationError(f"Acquisition id:{acquisition_id} not found")

    parser = OmeTiffParser(acquisition.location)
    acquisition_data = parser.get_acquisition_data()

    nuclei_channels = acquisition_data.get_image_stack_by_names(params.nuclei_channels)
    cytoplasm_channels = acquisition_data.get_image_stack_by_names(params.cytoplasm_channels)

    if params.preprocessing.channels_normalization == "minmax":
        nuclei_channels = normalize_by_minmax(nuclei_channels)
        cytoplasm_channels = normalize_by_minmax(cytoplasm_channels)
    elif params.preprocessing.channels_normalization == "zscore":
        nuclei_channels = normalize_by_zscore(nuclei_channels)
        cytoplasm_channels = normalize_by_zscore(cytoplasm_channels)

    nuclei_channels = np.nanmean(nuclei_channels, axis=0)
    cytoplasm_channels = np.nanmean(cytoplasm_channels, axis=0)

    # Combined together and expand to 4D
    im = np.stack((nuclei_channels, cytoplasm_channels), axis=-1)
    im = np.expand_dims(im, 0)

    app = Mesmer(model)
    segmentation_predictions = app.predict(
        im,
        batch_size=1,
        image_mpp=1.0,
        compartment=params.compartment,
        preprocess_kwargs=params.preprocessing.dict(),
        postprocess_kwargs_whole_cell=params.postprocessing.dict(),
        postprocess_kwargs_nuclear=params.postprocessing.dict(),
    )

    mask_filename = os.path.basename(acquisition.location).replace("ome.tiff", "mask.tiff")
    tifffile.imwrite(os.path.join(dataset.location, mask_filename), segmentation_predictions[0, :, :, 0])

    output = {
        "acquisition": acquisition,
        "mask_location": os.path.join(dataset.location, mask_filename),
        "object_numbers": None,
        "centroids_x": None,
        "centroids_y": None,
        "mean_intensities": {},
        "channel_names": acquisition_data.channel_names,
    }

    for c in acquisition_data.channel_names:
        d = measure.regionprops_table(
            label_image=segmentation_predictions[0, :, :, 0],
            intensity_image=acquisition_data.get_image_by_name(c),
            properties=("label", "centroid", "mean_intensity"),
        )
        if output["object_numbers"] is None:
            output["object_numbers"] = d.get("label")
            output["centroids_x"] = d.get("centroid-1")
            output["centroids_y"] = d.get("centroid-0")
        output["mean_intensities"][c] = d.get("mean_intensity")

    return output
