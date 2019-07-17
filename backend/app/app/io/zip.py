from __future__ import annotations

import csv
import logging
import os
import shutil
import zipfile
from pathlib import Path
from typing import Dict

import numpy as np
from imctools.io import mcdxmlparser
from imctools.io.imcacquisition import ImcAcquisition
from imctools.io.ometiffparser import OmetiffParser
from sqlalchemy.orm import Session

from app.core.errors import SlideImportError
from app.core.utils import timeit
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


@timeit
def import_zip(db: Session, uri: str, experiment_id: int):
    path = Path(uri)
    dir = path.parent
    with zipfile.ZipFile(path, 'r') as zip:
        zip.extractall(path.parent)

    subdirs = next(os.walk(dir))[1]
    if len(subdirs) > 0:
        dir = dir / subdirs[0]

    slide_meta_csv = _find_filename(dir, '_'.join([mcdxmlparser.SLIDE, 'meta.csv']))
    basename = slide_meta_csv.stem.replace(f'_{mcdxmlparser.SLIDE}_meta', '')

    slide_data = _load_meta_csv(slide_meta_csv)
    panorama_data = _load_meta_csv(dir / '_'.join([basename, mcdxmlparser.PANORAMA, 'meta.csv']))
    roi_data = _load_meta_csv(dir / '_'.join([basename, mcdxmlparser.ACQUISITIONROI, 'meta.csv']))
    roi_point_data = _load_meta_csv(dir / '_'.join([basename, mcdxmlparser.ROIPOINT, 'meta.csv']))
    acquisition_data = _load_meta_csv(dir / '_'.join([basename, mcdxmlparser.ACQUISITION, 'meta.csv']))
    channel_data = _load_meta_csv(dir / '_'.join([basename, mcdxmlparser.ACQUISITIONCHANNEL, 'meta.csv']))

    slide_map: Dict[str, Slide] = dict()
    for slide_item in slide_data.values():
        original_id = slide_item.get(mcdxmlparser.ID)
        metaname = f'{basename}_s{original_id}'
        item = slide_crud.get_by_metaname(db, experiment_id=experiment_id, metaname=metaname)
        if item:
            raise SlideImportError(
                f"The slide with metaname [{metaname}] already exists in the experiment [{experiment_id}]")

        xml_metadata_file: Path = dir / '_'.join([basename, "schema.xml"])
        with xml_metadata_file.open('rt') as f:
            original_metadata = f.read()

        slide = _import_slide(db, slide_item, original_metadata, experiment_id, basename)
        slide_map[str(slide.original_id)] = slide

        _copy_file(xml_metadata_file, slide.location)
        _copy_file(slide_meta_csv, slide.location)
        _copy_file(dir / '_'.join([basename, f's{slide.original_id}_slide.png']), slide.location)

    panorama_map: Dict[str, Panorama] = dict()
    for panorama_item in panorama_data.values():
        panorama = _import_panorama(db, panorama_item, slide_map.get(panorama_item.get(mcdxmlparser.SLIDEID)), basename)
        panorama_map[str(panorama.original_id)] = panorama
        _copy_file(dir / '_'.join([basename, f's{panorama.slide.original_id}_p{panorama.original_id}_pano.png']),
                   panorama.location)

    roi_map: Dict[str, ROI] = dict()
    for roi_item in roi_data.values():
        roi = _import_roi(db, roi_item, panorama_map.get(roi_item.get(mcdxmlparser.PANORAMAID)), basename)
        roi_map[str(roi.original_id)] = roi

    roi_point_map: Dict[str, ROIPoint] = dict()
    for roi_point_item in roi_point_data.values():
        roi_point = _import_roi_point(db, roi_point_item,
                                      roi_map.get(roi_point_item.get(mcdxmlparser.ACQUISITIONROIID)), basename)
        roi_point_map[str(roi_point.original_id)] = roi_point

    acquisition_map: Dict[str, Acquisition] = dict()
    for acquisition_item in acquisition_data.values():
        acquisition = _import_acquisition(db, acquisition_item,
                                          roi_map.get(acquisition_item.get(mcdxmlparser.ACQUISITIONROIID)), basename)
        acquisition_map[str(acquisition.original_id)] = acquisition

    for channel_item in channel_data.values():
        channel = _import_channel(db, channel_item, acquisition_map.get(channel_item.get(mcdxmlparser.ACQUISITIONID)),
                                  basename, dir)


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
    params = AcquisitionCreateModel(
        roi_id=roi.id,
        metaname=metaname,
        original_id=original_id,
        meta=meta,
    )
    acquisition = acquisition_crud.create(db, params=params)
    return acquisition


def _import_channel(db: Session, meta: Dict[str, str], acquisition: Acquisition, basename: str, dir: Path):
    original_id = meta.get(mcdxmlparser.ID)
    metaname = f'{basename}_s{acquisition.roi.panorama.slide.original_id}_p{acquisition.roi.panorama.original_id}_r{acquisition.roi.original_id}_a{acquisition.original_id}_c{original_id}'

    if meta.get(mcdxmlparser.CHANNELNAME) in ['X', 'Y', 'Z']:
        return

    metal = meta.get(mcdxmlparser.CHANNELNAME).replace('(', '').replace(')', '').strip()
    label = meta.get(mcdxmlparser.CHANNELLABEL).replace('(', '').replace(')', '').strip()
    mass = ''.join([m for m in metal if m.isdigit()])

    filename = dir / '_'.join([acquisition.metaname, 'ac.ome.tiff'])
    parser = OmetiffParser(filename)
    imc_acquisition: ImcAcquisition = parser.get_imc_acquisition()
    img = imc_acquisition.get_img_by_label(label)

    params = ChannelCreateModel(
        acquisition_id=acquisition.id,
        metaname=metaname,
        original_id=original_id,
        metal=metal,
        label=label,
        mass=mass,
        max_intensity=img.max(),
        min_intensity=img.min(),
        meta=meta
    )
    channel = channel_crud.create(db, params=params)
    np.save(os.path.join(channel.location, "origin.npy"), img)
    return channel


def _find_filename(dir: Path, ending: str) -> Path:
    filename = None
    for file in os.listdir(dir):
        if file.endswith(ending):
            filename = dir / file

    if not filename:
        raise SlideImportError(f"No file [{ending}] in zip folder [{dir}]")

    return filename


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


def _copy_file(src: Path, dst: str):
    return shutil.copy2(src, dst)
