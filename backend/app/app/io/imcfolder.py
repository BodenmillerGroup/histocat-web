from __future__ import annotations

import csv
import logging
import os
from pathlib import Path
from typing import Dict

from imctools.io import mcdxmlparser
from imctools.io.imcacquisition import ImcAcquisition
from imctools.io.ometiffparser import OmetiffParser
from sqlalchemy.orm import Session

from app.io.utils import copy_dir
from app.modules.acquisition import crud as acquisition_crud
from app.modules.acquisition.db import Acquisition
from app.modules.acquisition.models import AcquisitionCreateModel
from app.modules.channel import crud as channel_crud
from app.modules.channel.models import ChannelCreateModel
from app.modules.panorama import crud as panorama_crud
from app.modules.panorama.db import Panorama
from app.modules.panorama.models import PanoramaCreateModel
from app.modules.roi import crud as roi_crud
from app.modules.roi.db import ROI
from app.modules.roi.models import ROICreateModel
from app.modules.roi_point import crud as roi_point_crud
from app.modules.roi_point.db import ROIPoint
from app.modules.roi_point.models import ROIPointCreateModel
from app.modules.slide import crud as slide_crud
from app.modules.slide.db import Slide
from app.modules.slide.models import SlideCreateModel

logger = logging.getLogger(__name__)

SLIDE_META_CSV_ENDING = f"_{mcdxmlparser.SLIDE}{mcdxmlparser.META_CSV}"
PANORAMA_META_CSV_ENDING = f"_{mcdxmlparser.PANORAMA}{mcdxmlparser.META_CSV}"
ROI_META_CSV_ENDING = f"_{mcdxmlparser.ACQUISITIONROI}{mcdxmlparser.META_CSV}"
ROI_POINT_META_CSV_ENDING = f"_{mcdxmlparser.ROIPOINT}{mcdxmlparser.META_CSV}"
ACQUISITION_META_CSV_ENDING = f"_{mcdxmlparser.ACQUISITION}{mcdxmlparser.META_CSV}"
CHANNEL_META_CSV_ENDING = f"_{mcdxmlparser.ACQUISITIONCHANNEL}{mcdxmlparser.META_CSV}"


def import_imcfolder(db: Session, schema_filename: str, experiment_id: int):
    """
    Import slides from the folder compatible with IMC pipeline
    """

    schema_path = Path(schema_filename)
    src_folder = schema_path.parent
    basename = schema_path.stem.replace('_schema', '')

    slide_csv_filename = src_folder / f"{basename}{SLIDE_META_CSV_ENDING}"
    slide_data = _load_meta_csv(slide_csv_filename)

    panorama_csv_filename = src_folder / f"{basename}{PANORAMA_META_CSV_ENDING}"
    panorama_data = _load_meta_csv(panorama_csv_filename)

    roi_csv_filename = src_folder / f"{basename}{ROI_META_CSV_ENDING}"
    roi_data = _load_meta_csv(roi_csv_filename)

    roi_point_csv_filename = src_folder / f"{basename}{ROI_POINT_META_CSV_ENDING}"
    roi_point_data = _load_meta_csv(roi_point_csv_filename)

    acquisition_csv_filename = src_folder / f"{basename}{ACQUISITION_META_CSV_ENDING}"
    acquisition_data = _load_meta_csv(acquisition_csv_filename)

    channel_csv_filename = src_folder / f"{basename}{CHANNEL_META_CSV_ENDING}"
    channel_data = _load_meta_csv(channel_csv_filename)

    slide_map: Dict[str, Slide] = dict()
    for slide_item in slide_data.values():
        original_id = slide_item.get(mcdxmlparser.ID)
        metaname = f'{basename}_s{original_id}'
        item = slide_crud.get_by_metaname(db, experiment_id=experiment_id, metaname=metaname)
        if item:
            logger.warn(f"The slide with metaname [{metaname}] already exists in the experiment [{experiment_id}]")
            return

        with open(schema_filename, 'rt') as f:
            original_metadata = f.read()

        slide = _import_slide(db, slide_item, original_metadata, experiment_id, basename)
        slide_map[str(slide.original_id)] = slide

        origin_location = os.path.join(slide.location, 'origin')
        copy_dir(src_folder, origin_location)

    panorama_map: Dict[str, Panorama] = dict()
    for panorama_item in panorama_data.values():
        slide = slide_map.get(panorama_item.get(mcdxmlparser.SLIDEID))
        panorama = _import_panorama(db, panorama_item, slide, basename)
        panorama_map[str(panorama.original_id)] = panorama

    roi_map: Dict[str, ROI] = dict()
    for roi_item in roi_data.values():
        panorama = panorama_map.get(roi_item.get(mcdxmlparser.PANORAMAID))
        roi = _import_roi(db, roi_item, panorama, basename)
        roi_map[str(roi.original_id)] = roi

    roi_point_map: Dict[str, ROIPoint] = dict()
    for roi_point_item in roi_point_data.values():
        roi = roi_map.get(roi_point_item.get(mcdxmlparser.ACQUISITIONROIID))
        roi_point = _import_roi_point(db, roi_point_item, roi, basename)
        roi_point_map[str(roi_point.original_id)] = roi_point

    acquisition_map: Dict[str, Acquisition] = dict()
    for acquisition_item in acquisition_data.values():
        roi = roi_map.get(acquisition_item.get(mcdxmlparser.ACQUISITIONROIID))
        acquisition = _import_acquisition(db, acquisition_item, roi, basename)
        acquisition_map[str(acquisition.original_id)] = acquisition

    for channel_item in channel_data.values():
        acquisition = acquisition_map.get(channel_item.get(mcdxmlparser.ACQUISITIONID))
        _import_channel(db, channel_item, acquisition, basename)


def _import_slide(db: Session, meta: Dict[str, str], original_metadata: str, experiment_id: int, basename: str):
    original_id = meta.get(mcdxmlparser.ID)
    metaname = f'{basename}_s{original_id}'
    params = SlideCreateModel(
        experiment_id=experiment_id,
        metaname=metaname,
        original_id=original_id,
        original_metadata=original_metadata,
        meta=meta,
    )
    slide = slide_crud.create(db, params=params)
    return slide


def _import_panorama(db: Session, meta: Dict[str, str], slide: Slide, basename: str):
    original_id = meta.get(mcdxmlparser.ID)
    metaname = f'{basename}_s{slide.original_id}_p{original_id}'
    params = PanoramaCreateModel(
        slide_id=slide.id,
        metaname=metaname,
        original_id=original_id,
        meta=meta,
    )
    panorama = panorama_crud.create(db, params=params)
    return panorama


def _import_roi(db: Session, meta: Dict[str, str], panorama: Panorama, basename: str):
    original_id = meta.get(mcdxmlparser.ID)
    metaname = f'{basename}_s{panorama.slide.original_id}_p{panorama.original_id}_r{original_id}'
    params = ROICreateModel(
        panorama_id=panorama.id,
        metaname=metaname,
        original_id=original_id,
        meta=meta,
    )
    roi = roi_crud.create(db, params=params)
    return roi


def _import_roi_point(db: Session, meta: Dict[str, str], roi: ROI, basename: str):
    original_id = meta.get(mcdxmlparser.ID)
    metaname = f'{basename}_s{roi.panorama.slide.original_id}_p{roi.panorama.original_id}_r{roi.original_id}_rp{original_id}'
    params = ROIPointCreateModel(
        roi_id=roi.id,
        metaname=metaname,
        original_id=original_id,
        meta=meta,
    )
    roi_point = roi_point_crud.create(db, params=params)
    return roi_point


def _import_acquisition(db: Session, meta: Dict[str, str], roi: ROI, basename: str):
    original_id = meta.get(mcdxmlparser.ID)
    metaname = f'{basename}_s{roi.panorama.slide.original_id}_p{roi.panorama.original_id}_r{roi.original_id}_a{original_id}'
    origin_location = os.path.join(roi.panorama.slide.location, 'origin')
    location = os.path.join(origin_location, f'{metaname}_ac.ome.tiff')
    params = AcquisitionCreateModel(
        roi_id=roi.id,
        metaname=metaname,
        original_id=original_id,
        location=location,
        meta=meta,
    )
    acquisition = acquisition_crud.create(db, params=params)
    return acquisition


def _import_channel(db: Session, meta: Dict[str, str], acquisition: Acquisition, basename: str):
    original_id = meta.get(mcdxmlparser.ID)
    metaname = f'{basename}_s{acquisition.roi.panorama.slide.original_id}_p{acquisition.roi.panorama.original_id}_r{acquisition.roi.original_id}_a{acquisition.original_id}_c{original_id}'

    if meta.get(mcdxmlparser.CHANNELNAME) in ['X', 'Y', 'Z']:
        return

    metal = meta.get(mcdxmlparser.CHANNELNAME).replace('(', '').replace(')', '').strip()
    label = meta.get(mcdxmlparser.CHANNELLABEL).replace('(', '').replace(')', '').strip()
    mass = ''.join([m for m in metal if m.isdigit()])

    parser = OmetiffParser(acquisition.location)
    imc_acquisition: ImcAcquisition = parser.get_imc_acquisition()
    img = imc_acquisition.get_img_by_label(label)

    params = ChannelCreateModel(
        acquisition_id=acquisition.id,
        metaname=metaname,
        original_id=original_id,
        metal=metal,
        label=label,
        mass=mass,
        location=acquisition.location,
        max_intensity=img.max(),
        min_intensity=img.min(),
        meta=meta
    )
    channel = channel_crud.create(db, params=params)
    return channel


def _load_meta_csv(filepath: Path) -> Dict[str, Dict[str, str]]:
    with filepath.open('rt') as f:
        reader = csv.reader(f)
        line_count = 0
        header = None
        data = dict()
        for row in reader:
            if line_count == 0:
                header = row
            else:
                meta = dict()
                for i, v in enumerate(row):
                    meta[header[i]] = v
                id = meta.get(mcdxmlparser.ID)
                data[id] = meta
            line_count += 1
        return data
