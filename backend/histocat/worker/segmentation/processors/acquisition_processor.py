import os

import cv2
import numpy as np
import tifffile
from histocat.core.errors import SegmentationError
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from skimage import exposure, measure
from skimage.filters import gaussian
from skimage.transform import rescale
from sqlalchemy.orm import Session
from deepcell.datasets import multiplex_tissue
from deepcell.utils.plot_utils import create_rgb_image
from deepcell.applications import MultiplexSegmentation
from deepcell.utils.plot_utils import make_outline_overlay
from imctoolkit import MultichannelImage, ImageSingleCellData

from histocat.core.acquisition import service as acquisition_service
from histocat.core.segmentation.dto import SegmentationSubmissionDto
from histocat.core.utils import timeit


@timeit
def process_acquisition(db: Session, acquisition_id: int, params: SegmentationSubmissionDto, model):
    app = MultiplexSegmentation(model)
    print('Training Resolution:', app.model_mpp, 'microns per pixel')

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

    # rescale intensities nucleus
    # for j in range(0, 1):
    #     p0, p1 = np.percentile(im[:, :, j], (0, params.upper_limit))
    #     im[:, :, j] = exposure.rescale_intensity(im[:, :, j], in_range=(p0, p1))



    print(im.shape)
    segmentation_predictions = app.predict(im, image_mpp=0.5, compartment="whole-cell")
    rgb_images = create_rgb_image(im, channel_colors=['green', 'blue'])
    overlay_data = make_outline_overlay(rgb_data=rgb_images, predictions=segmentation_predictions)
    print(overlay_data.shape)

    filename = os.path.basename(acquisition.location).replace("ome.tiff", "overlay.tiff")
    tifffile.imwrite(os.path.join("/data/inbox", filename), overlay_data[0, :, :, :])

    filename = os.path.basename(acquisition.location).replace("ome.tiff", "mask.tiff")
    tifffile.imwrite(os.path.join("/data/inbox", filename), segmentation_predictions[0, :, :, 0])

    output = {}
    for c in acquisition_data.channel_names:
        d = measure.regionprops_table(label_image=segmentation_predictions[0, :, :, 0], intensity_image=acquisition_data.get_image_by_name(c), properties=('label', 'centroid', 'mean_intensity'))
        output[c] = d

    # img = MultichannelImage.read_tiff(acquisition.location, ome_channel_name_attr="Fluor")
    # d = ImageSingleCellData(img, segmentation_predictions[0, :, :, 0])
    # c = d.cell_centroids
    # m = d.mean_intensities
    pass
