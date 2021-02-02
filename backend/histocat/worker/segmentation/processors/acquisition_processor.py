import os

import numpy as np
import tifffile
from deepcell.applications import MultiplexSegmentation
from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from skimage import io, measure
from sqlalchemy.orm import Session

from histocat.core.acquisition import service as acquisition_service
from histocat.core.dataset.models import DatasetModel
from histocat.core.errors import SegmentationError
from histocat.core.segmentation.dto import SegmentationSubmissionDto
from histocat.core.utils import timeit


@timeit
def process_acquisition(
    db: Session, acquisition_id: int, params: SegmentationSubmissionDto, model, dataset: DatasetModel
):
    app = MultiplexSegmentation(model)

    acquisition = acquisition_service.get_by_id(db, acquisition_id)
    if not acquisition:
        raise SegmentationError(f"Acquisition id:{acquisition_id} not found")

    parser = OmeTiffParser(acquisition.location)
    acquisition_data = parser.get_acquisition_data()

    IA_stack = np.zeros(
        (acquisition_data.image_data.shape[1], acquisition_data.image_data.data.shape[2], 2)
    )  # Channels last

    # sum up nuclear markers
    for c in params.nuclei_channels:
        img = acquisition_data.get_image_by_name(c)
        IA_stack[:, :, 0] = IA_stack[:, :, 0] + img

    # sum up cytoplasmic/membranous markers
    for c in params.cytoplasm_channels:
        img = acquisition_data.get_image_by_name(c)
        IA_stack[:, :, 1] = IA_stack[:, :, 1] + img

    im = np.stack((IA_stack[:, :, 0], IA_stack[:, :, 1]), axis=-1)
    im = np.expand_dims(im, 0)

    segmentation_predictions = app.predict(
        im,
        batch_size=1,
        image_mpp=1.0,
        compartment="whole-cell",
        preprocess_kwargs=params.preprocessing.dict(),
        postprocess_kwargs_whole_cell=params.postprocessing.dict(),
    )
    rgb_images = create_rgb_image(im, channel_colors=["green", "blue"])
    overlay_data = make_outline_overlay(rgb_data=rgb_images, predictions=segmentation_predictions)

    filename = os.path.basename(acquisition.location).replace("ome.tiff", "overlay.png")
    io.imsave(os.path.join(dataset.location, filename), overlay_data[0, :, :, :])

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
