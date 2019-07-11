from __future__ import annotations

import os
from xml.etree import ElementTree

import numpy as np
from imctools.io.mcdparser import McdParser
from sqlalchemy.orm import Session

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


def _import_slide(db: Session, item, original_metadata: str, experiment_id: int):
    params = SlideCreateModel(
        experiment_id=experiment_id,
        metaname=item.metaname,
        original_id=item.id,
        original_metadata=original_metadata,
        meta=item.properties,
    )
    slide = slide_crud.create(db, params=params)
    return slide


def _import_panorama(db: Session, item, slide_id: int):
    params = PanoramaCreateModel(
        slide_id=slide_id,
        metaname=item.metaname,
        original_id=item.id,
        meta=item.properties,
    )
    panorama = panorama_crud.create(db, params=params)
    return panorama


def _import_roi(db: Session, item, panorama_id: int):
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


@timeit
def import_mcd(db: Session, uri: str, experiment_id: int):
    with McdParser(uri) as mcd:
        slide_item = mcd.meta.objects["Slide"]["0"]
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
