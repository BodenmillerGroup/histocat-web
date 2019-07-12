import csv
import logging
import os
import re
import shutil
from collections import namedtuple
from typing import Dict

import numpy as np
from imctools.io.tiffwriter import TiffWriter
from sqlalchemy.orm import Session

from app.core.errors import DatasetInputError
from app.modules.acquisition import crud as acquisition_crud
from app.modules.dataset.db import Dataset
from app.modules.panorama.db import Panorama
from app.modules.roi.db import ROI
from app.modules.slide.db import Slide

logger = logging.getLogger(__name__)

FileLink = namedtuple('FileLink', ['src', 'dst'])


def prepare_dataset(db: Session, dataset: Dataset):
    # output folder where the output files will be saved
    folder_base = dataset.location

    # parameters for resizing the images for ilastik
    folder_analysis = os.path.join(folder_base, 'tiffs')
    folder_ilastik = os.path.join(folder_base, 'ilastik')
    folder_ome = os.path.join(folder_base, 'ometiff')
    folder_cp = os.path.join(folder_base, 'cpout')
    folder_histocat = os.path.join(folder_base, 'histocat')
    folder_uncertainty = os.path.join(folder_base, 'uncertainty')

    for folder in [folder_base, folder_analysis, folder_ilastik, folder_ome, folder_cp, folder_histocat,
                   folder_uncertainty]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    input = dataset.input
    acquisition_ids = input.get('acquisition_ids')
    if not acquisition_ids or len(acquisition_ids) == 0:
        raise DatasetInputError(f'Dataset [{dataset.id}]: Acquisitions are not selected')
    metals = input.get('metals')
    if not metals or len(metals) == 0:
        raise DatasetInputError(f'Dataset [{dataset.id}]: Metals are not selected')
    channel_settings = input.get('channel_settings')

    slide_csv: Dict[int, dict] = dict()
    panorama_csv: Dict[int, dict] = dict()
    roi_csv: Dict[int, dict] = dict()
    acquisition_csv: Dict[int, dict] = dict()
    channel_csv: Dict[int, dict] = dict()
    roi_point_csv: Dict[int, dict] = dict()

    short_slide_name = None
    file_links = []

    for acquisition_id in acquisition_ids:
        acquisition = acquisition_crud.get(db, id=acquisition_id)
        roi: ROI = acquisition.roi
        panorama: Panorama = roi.panorama
        slide: Slide = panorama.slide

        short_slide_name = re.sub('_s\d+$', '', slide.metaname)

        slide_csv[slide.id] = slide.meta
        panorama_csv[panorama.id] = panorama.meta
        roi_csv[roi.id] = roi.meta
        acquisition_csv[acquisition.id] = acquisition.meta

        for roi_point in roi.roi_points:
            roi_point_csv[roi_point.id] = roi_point.meta

        channel_arrays = []
        channel_names = []
        channel_fluors = []
        for channel in acquisition.channels:
            if channel.metal in metals:
                channel_csv[channel.id] = channel.meta
                channel_arrays.append(np.load(os.path.join(channel.location, "origin.npy")))
                channel_names.append(channel.label)
                channel_fluors.append(channel.metal)

        img_stack = np.stack(channel_arrays).swapaxes(2, 0)

        tiff_writer = TiffWriter(
            os.path.join(folder_ome, acquisition.metaname + '_ac.ome.tiff'),
            img_stack=img_stack,
            channel_name=channel_names,
            original_description=slide.original_metadata,
            fluor=channel_fluors,
        )

        tiff_writer.save_image(
            mode='ome',
            compression=0,
            bigtiff=True,
        )

        file_links.extend([
            FileLink(
                src=os.path.join(slide.location, short_slide_name + '_schema.xml'),
                dst=os.path.join(folder_ome, short_slide_name + '_schema.xml')
            ),
            FileLink(
                src=os.path.join(slide.location, slide.metaname + '_slide.png'),
                dst=os.path.join(folder_ome, slide.metaname + '_slide.png')
            ),
            FileLink(
                src=os.path.join(panorama.location, panorama.metaname + '_pano.png'),
                dst=os.path.join(folder_ome, panorama.metaname + '_pano.png')
            ),
        ])

    for link in file_links:
        if not os.path.exists(link.dst):
            shutil.copy2(link.src, link.dst)

    _save_meta_csv(slide_csv, os.path.join(folder_ome, short_slide_name + '_Slide_meta.csv'))
    _save_meta_csv(panorama_csv, os.path.join(folder_ome, short_slide_name + '_Panorama_meta.csv'))
    _save_meta_csv(roi_csv, os.path.join(folder_ome, short_slide_name + '_AcquisitionROI_meta.csv'))
    _save_meta_csv(acquisition_csv, os.path.join(folder_ome, short_slide_name + '_Acquisition_meta.csv'))
    _save_meta_csv(channel_csv, os.path.join(folder_ome, short_slide_name + '_AcquisitionChannel_meta.csv'))
    _save_meta_csv(roi_point_csv, os.path.join(folder_ome, short_slide_name + '_ROIPoint_meta.csv'))


def _save_meta_csv(items: Dict[int, dict], filename: str):
    """
    Writes the xml data as csv tables

    """
    with open(filename, 'w') as f:
        cols = next(iter(items.values())).keys()
        writer = csv.DictWriter(f, cols)
        writer.writeheader()
        for row in items.values():
            writer.writerow(row)
