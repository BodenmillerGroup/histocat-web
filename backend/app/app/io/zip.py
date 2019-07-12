from __future__ import annotations

import csv
import os
from typing import Dict
from xml.etree import ElementTree
import zipfile
import shutil
from pathlib import Path

import numpy as np
from sqlalchemy.orm import Session

from app.core.errors import SlideImportError
from app.core.utils import timeit
from app.modules.acquisition import crud as acquisition_crud
from app.modules.acquisition.models import AcquisitionCreateModel
from app.modules.channel import crud as channel_crud
from app.modules.channel.models import ChannelCreateModel
from app.modules.panorama import crud as panorama_crud
from app.modules.panorama.models import PanoramaCreateModel
from app.modules.roi import crud as roi_crud
from app.modules.roi.models import ROICreateModel
from app.modules.roi_point import crud as roi_point_crud
from app.modules.roi_point.models import ROIPointCreateModel
from app.modules.slide import crud as slide_crud
from app.modules.slide.models import SlideCreateModel


@timeit
def import_zip(db: Session, uri: str, experiment_id: int):
    path = Path(uri)
    dir = path.parent
    basename = path.stem
    with zipfile.ZipFile(path, 'r') as zip:
        zip.extractall(path.parent)

    subdirs = next(os.walk(dir))[1]
    if len(subdirs) > 0:
        dir = dir / subdirs[0]

    slide_data = _load_meta_csv(dir, "_Slide_meta.csv")
    panorama_data = _load_meta_csv(dir, "_Panorama_meta.csv")
    roi_data = _load_meta_csv(dir, "_AcquisitionROI_meta.csv")
    roi_point_data = _load_meta_csv(dir, "_ROIPoint_meta.csv")
    acquisition_data = _load_meta_csv(dir, "_Acquisition_meta.csv")
    channel_data = _load_meta_csv(dir, "_AcquisitionChannel_meta.csv")

    for slide_item in slide_data.values():
        uid = slide_item.get('UID')
        item = slide_crud.get_by_uid(db, uid=uid)
        if item:
            raise SlideImportError(f"The slide with UID [{uid}] already exists in the experiment [{experiment_id}]")

        xml_metadata_filename = _find_filename(dir, "_schema.xml")
        original_metadata = None
        with open(xml_metadata_filename, 'rt') as f:
            original_metadata = f.read()

        slide = _import_slide(db, slide_item, original_metadata, experiment_id, basename)
        shutil.copy2(xml_metadata_filename, os)



    with McdParser(uri) as mcd:
        slide_item = mcd.meta.objects["Slide"]["0"]
        uid = slide_item.properties.get('UID')
        item = slide_crud.get_by_uid(db, uid=uid)
        if item:
            raise SlideImportError(f"The slide with UID [{uid}] already exists in the experiment [{experiment_id}]")
        original_metadata = ElementTree.tostring(mcd.xml, encoding="utf8", method="xml")
        slide = _import_slide(db, slide_item, original_metadata, experiment_id)
        mcd.meta.save_meta_xml(slide.location)
        mcd.meta.save_meta_csv(slide.location)
        mcd.save_slideimage(slide_item.id, slide.location)

        for panorama_item in slide_item.childs["Panorama"].values():
            panorama = _import_panorama(db, panorama_item, slide.id)
            mcd.save_panorama(panorama_item.id, panorama.location)

            if "AcquisitionROI" in panorama_item.childs:
                for roi_item in panorama_item.childs["AcquisitionROI"].values():
                    roi = _import_roi(db, roi_item, panorama.id)

                    for roi_point_item in roi_item.childs["ROIPoint"].values():
                        roi_point = _import_roi_point(db, roi_point_item, roi.id)

                    for acquisition_item in roi_item.childs["Acquisition"].values():
                        acquisition = _import_acquisition(db, acquisition_item, roi.id)
                        mcd.save_acquisition_bfimage_before(acquisition_item.id, acquisition.location)
                        mcd.save_acquisition_bfimage_after(acquisition_item.id, acquisition.location)

                        imc_acquisition = mcd.get_imc_acquisition(acquisition_item.id)
                        for channel_item in acquisition_item.childs["AcquisitionChannel"].values():
                            channel = _import_channel(db, channel_item, imc_acquisition, acquisition.id)


def _import_slide(db: Session, meta, original_metadata: str, experiment_id: int, basename: str):
    original_id = meta.get('ID')
    metaname = f'{basename}_s{original_id}'
    params = SlideCreateModel(
        experiment_id=experiment_id,
        metaname=metaname,
        original_id=original_id,
        uid=meta.get('UID'),
        original_metadata=original_metadata,
        meta=meta,
    )
    slide = slide_crud.create(db, params=params)
    return slide


def _import_panorama(db: Session, item, slide_id: int, basename: str):
    params = PanoramaCreateModel(
        slide_id=slide_id,
        metaname=item.metaname,
        original_id=item.id,
        meta=item.properties,
    )
    panorama = panorama_crud.create(db, params=params)
    return panorama


def _import_roi(db: Session, item, panorama_id: int, basename: str):
    params = ROICreateModel(
        panorama_id=panorama_id,
        metaname=item.metaname,
        original_id=item.id,
        meta=item.properties,
    )
    roi = roi_crud.create(db, params=params)
    return roi


def _import_roi_point(db: Session, item, roi_id: int):
    params = ROIPointCreateModel(
        roi_id=roi_id,
        metaname=item.metaname,
        original_id=item.id,
        meta=item.properties,
    )
    roi_point = roi_point_crud.create(db, params=params)
    return roi_point


def _import_acquisition(db: Session, item, roi_id: int):
    params = AcquisitionCreateModel(
        roi_id=roi_id,
        metaname=item.metaname,
        original_id=item.id,
        meta=item.properties,
    )
    acquisition = acquisition_crud.create(db, params=params)
    return acquisition


def _import_channel(db: Session, item, imc_acquisition, acquisition_id: int):
    if item.properties['ChannelName'] in ['X', 'Y', 'Z']:
        return

    metal = item.properties['ChannelName'].replace('(', '').replace(')', '').strip()
    label = item.properties['ChannelLabel'].replace('(', '').replace(')', '').strip()
    mass = ''.join([m for m in metal if m.isdigit()])
    img = imc_acquisition.get_img_by_label(label)

    params = ChannelCreateModel(
        acquisition_id=acquisition_id,
        metaname=item.metaname,
        original_id=item.id,
        metal=metal,
        label=label,
        mass=mass,
        max_intensity=img.max(),
        min_intensity=img.min(),
        meta=item.properties
    )
    channel = channel_crud.create(db, params=params)
    np.save(os.path.join(channel.location, "origin.npy"), img)
    return channel


def _find_filename(dir: str, ending: str) -> str:
    filename = None
    for file in os.listdir(dir):
        if file.endswith(ending):
            filename = os.path.join(dir, file)

    if not filename:
        raise SlideImportError(f"No file [{ending}] in zip folder [{dir}]")

    return filename


def _load_meta_csv(dir: str, ending: str) -> Dict[str, Dict[str, str]]:
    filename = _find_filename(dir, ending)
    with open(filename, 'rt') as f:
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
                id = meta.get('ID')
                data[id] = meta
            line_count += 1
        return data


def _copy_file(src_dir: str, dst_dir: str, ending: str):
    src = _find_filename(src_dir, ending)
    filename = os.path.basename()
    shutil.copy2(src)
