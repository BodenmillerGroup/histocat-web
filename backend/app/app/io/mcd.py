from __future__ import annotations

import os

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


def _import_slide(db: Session, item, experiment_id: int):
    params = SlideCreateModel(
        experiment_id=experiment_id,
        metaname=item.metaname,
        original_id=item.id,
        uid=item.properties["UID"],
        description=item.properties["Description"],
        filename=item.properties["Filename"],
        slide_type=item.properties["SlideType"],
        width_um=item.properties["WidthUm"],
        height_um=item.properties["HeightUm"],
        image_end_offset=item.properties["ImageEndOffset"],
        image_start_offset=item.properties["ImageStartOffset"],
        image_file=item.properties["ImageFile"],
        meta=item.properties,
    )
    slide = slide_crud.create(db, params=params)
    return slide


def _import_panorama(db: Session, item, slide_id: int):
    params = PanoramaCreateModel(
        slide_id=slide_id,
        metaname=item.metaname,
        original_id=item.id,
        description=item.properties["Description"],
        slide_x1_pos_um=item.properties["SlideX1PosUm"],
        slide_y1_pos_um=item.properties["SlideY1PosUm"],
        slide_x2_pos_um=item.properties["SlideX2PosUm"],
        slide_y2_pos_um=item.properties["SlideY2PosUm"],
        slide_x3_pos_um=item.properties["SlideX3PosUm"],
        slide_y3_pos_um=item.properties["SlideY3PosUm"],
        slide_x4_pos_um=item.properties["SlideX4PosUm"],
        slide_y4_pos_um=item.properties["SlideY4PosUm"],
        image_end_offset=item.properties["ImageEndOffset"],
        image_start_offset=item.properties["ImageStartOffset"],
        pixel_width=item.properties["PixelWidth"],
        pixel_height=item.properties["PixelHeight"],
        image_format=item.properties["ImageFormat"],
        pixel_scale_coef=item.properties["PixelScaleCoef"],
        meta=item.properties,
    )
    panorama = panorama_crud.create(db, params=params)
    return panorama


def _import_roi(db: Session, item, panorama_id: int):
    params = ROICreateModel(
        panorama_id=panorama_id,
        metaname=item.metaname,
        original_id=item.id,
        roi_type=item.properties["ROIType"],
    )
    roi = roi_crud.create(db, params=params)
    return roi


def _import_roi_point(db: Session, item, roi_id: int):
    params = ROIPointCreateModel(
        roi_id=roi_id,
        metaname=item.metaname,
        original_id=item.id,
        order_number=item.properties["OrderNumber"],
        slide_x_pos_um=item.properties["SlideXPosUm"],
        slide_y_pos_um=item.properties["SlideYPosUm"],
        panorama_pixel_x_pos=item.properties["PanoramaPixelXPos"],
        panorama_pixel_y_pos=item.properties["PanoramaPixelYPos"],
    )
    roi_point = roi_point_crud.create(db, params=params)
    return roi_point


def _import_acquisition(db: Session, item, roi_id: int):
    params = AcquisitionCreateModel(
        roi_id=roi_id,
        metaname=item.metaname,
        original_id=item.id,
        description=item.properties["Description"],
        order_number=item.properties["OrderNumber"],
        ablation_power=item.properties["AblationPower"],
        ablation_distance_between_shots_x=item.properties["AblationDistanceBetweenShotsX"],
        ablation_distance_between_shots_y=item.properties["AblationDistanceBetweenShotsY"],
        ablation_frequency=item.properties["AblationFrequency"],
        signal_type=item.properties["SignalType"],
        dual_count_start=item.properties["DualCountStart"],
        data_start_offset=item.properties["DataStartOffset"],
        data_end_offset=item.properties["DataEndOffset"],
        start_timestamp=item.properties["StartTimeStamp"],
        end_timestamp=item.properties["EndTimeStamp"],
        after_ablation_image_end_offset=item.properties["AfterAblationImageEndOffset"],
        after_ablation_image_start_offset=item.properties["AfterAblationImageStartOffset"],
        before_ablation_image_end_offset=item.properties["BeforeAblationImageEndOffset"],
        before_ablation_image_start_offset=item.properties["BeforeAblationImageStartOffset"],
        roi_start_x_pos_um=item.properties["ROIStartXPosUm"],
        roi_start_y_pos_um=item.properties["ROIStartYPosUm"],
        roi_end_x_pos_um=item.properties["ROIEndXPosUm"],
        roi_end_y_pos_um=item.properties["ROIEndYPosUm"],
        movement_type=item.properties["MovementType"],
        segment_data_format=item.properties["SegmentDataFormat"],
        value_bytes=item.properties["ValueBytes"],
        max_y=item.properties["MaxY"],
        max_x=item.properties["MaxX"],
        plume_start=item.properties["PlumeStart"],
        plume_end=item.properties["PlumeEnd"],
        template=item.properties["Template"],
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
        slide = _import_slide(db, slide_item, experiment_id)
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
