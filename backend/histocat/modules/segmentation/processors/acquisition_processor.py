import os

import cv2
import numpy as np
import tifffile
from fastapi import HTTPException
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from skimage import exposure
from skimage.filters import gaussian
from skimage.transform import rescale, resize
from sqlalchemy.orm import Session
from starlette import status
from tensorflow import keras
import cellprofiler_core.image
import cellprofiler_core.object
import cellprofiler_core.pipeline
import cellprofiler_core.preferences
import cellprofiler_core.workspace
from cellprofiler_core import preferences as cp_preferences
from cellprofiler_core import pipeline as cp_pipeline

from histocat.core.utils import timeit
from histocat.modules.acquisition import service as acquisition_service
from histocat.modules.model.models import ModelModel
from histocat.modules.segmentation.dto import SegmentationSubmissionDto

cp_preferences.set_headless()


@timeit
def process_acquisition(db: Session, acquisition_id: int, params: SegmentationSubmissionDto, model):
    acquisition = acquisition_service.get_by_id(db, acquisition_id)
    if not acquisition:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Acquisition id:{acquisition_id} not found")

    parser = OmeTiffParser(acquisition.location)
    acquisition_data = parser.get_acquisition_data()

    IA_stack = np.zeros(
        (acquisition_data.image_data.shape[1], acquisition_data.image_data.data.shape[2], 3)
    )  # Channels last

    # sum up nuclear markers
    for c in params.nuclei_channels:
        img = acquisition_data.get_image_by_name(c)
        IA_stack[:, :, 0] = IA_stack[:, :, 0] + img

    # sum up cytoplasmic/membranous markers
    for c in params.cytoplasm_channels:
        img = acquisition_data.get_image_by_name(c)
        IA_stack[:, :, 1] = IA_stack[:, :, 1] + img

    # resize by the given upscaling factor the RBG image
    IA_stack = rescale(IA_stack, params.scaling_factor, multichannel=True, anti_aliasing=False, mode="reflect")
    IA_stack = np.array(IA_stack, dtype="uint16") / 65535

    # rescale intensities nucleus
    for j in range(0, 2):
        p0, p1 = np.percentile(IA_stack[:, :, j], (0, params.upper_limit))
        IA_stack[:, :, j] = exposure.rescale_intensity(IA_stack[:, :, j], in_range=(p0, p1))

    # generate probabilities
    tiles_size = model.input.shape[2]
    channels = model.input.shape[3]

    filename = os.path.basename(acquisition.location).replace("ome.tiff", "ome.Probabilities.tiff")

    # crop the image
    imgwidth, imgheight = IA_stack.shape[0], IA_stack.shape[1]
    fovs_x = int(np.ceil(imgwidth / tiles_size))
    fovs_y = int(np.ceil(imgheight / tiles_size))

    pred = np.zeros((fovs_x * tiles_size, fovs_y * tiles_size, channels))  # depends on the model

    for i in range(
        0, imgheight, int(tiles_size * (1 - params.overlap))
    ):  # for the moment overlap the last bit of tile size
        for j in range(0, imgwidth, int(tiles_size * (1 - params.overlap))):
            frame = IA_stack[j : j + tiles_size, i : i + tiles_size, :]
            # frame = np.array(frame,dtype="uint16")
            xd = tiles_size - frame.shape[0]
            yd = tiles_size - frame.shape[1]
            tile = cv2.copyMakeBorder(frame, 0, xd, 0, yd, cv2.BORDER_CONSTANT)  # padding for last tiles
            # predict classes
            prediction = model.predict(np.expand_dims(tile, axis=0), verbose=0)
            prediction = prediction[0, :, :, :]
            x = pred[j : j + tiles_size, i : i + tiles_size, :].shape[0]
            y = pred[j : j + tiles_size, i : i + tiles_size, :].shape[1]
            pred[j : j + tiles_size, i : i + tiles_size, :] = prediction[
                0:x, 0:y,
            ]  # cropped prediction

    # crop prediction to original image size
    pred = pred[0:imgwidth, 0:imgheight, :]
    pred = gaussian(pred, sigma=params.filter_settings.sigma, multichannel=True)  ##minor gaussian blur

    # convert to uint16
    # rearrange the mask into RGB
    pred = np.array(pred * 65536, dtype="uint16")  # 16 bits for CP

    tifffile.imwrite(os.path.join("/data/inbox", filename), np.array(pred, dtype="uint16"))

    # Run CellProfiler pipeline
    pipeline = cp_pipeline.Pipeline()
    pipeline.load("/data/inbox/2_segment_ilastik.cppipe")
