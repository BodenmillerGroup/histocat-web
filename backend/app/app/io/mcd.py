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


def _import_slide(db: Session, meta: dict, experiment_id: int):
    params = SlideCreateModel(
        experiment_id=experiment_id,
        uid=meta["UID"],
        description=meta["Description"],
        filename=meta["Filename"],
        slide_type=meta["SlideType"],
        width_um=meta["WidthUm"],
        height_um=meta["HeightUm"],
        image_end_offset=meta["ImageEndOffset"],
        image_start_offset=meta["ImageStartOffset"],
        image_file=meta["ImageFile"],
        meta=meta,
    )
    slide = slide_crud.create(db, params=params)
    return slide


def _import_panorama(db: Session, meta: dict, slide_id: int):
    params = PanoramaCreateModel(
        slide_id=slide_id,
        description=meta["Description"],
        slide_x1_pos_um=meta["SlideX1PosUm"],
        slide_y1_pos_um=meta["SlideY1PosUm"],
        slide_x2_pos_um=meta["SlideX2PosUm"],
        slide_y2_pos_um=meta["SlideY2PosUm"],
        slide_x3_pos_um=meta["SlideX3PosUm"],
        slide_y3_pos_um=meta["SlideY3PosUm"],
        slide_x4_pos_um=meta["SlideX4PosUm"],
        slide_y4_pos_um=meta["SlideY4PosUm"],
        image_end_offset=meta["ImageEndOffset"],
        image_start_offset=meta["ImageStartOffset"],
        pixel_width=meta["PixelWidth"],
        pixel_height=meta["PixelHeight"],
        image_format=meta["ImageFormat"],
        pixel_scale_coef=meta["PixelScaleCoef"],
        meta=meta,
    )
    panorama = panorama_crud.create(db, params=params)
    return panorama


def _import_roi(db: Session, meta: dict, panorama_id: int):
    params = ROICreateModel(
        panorama_id=panorama_id,
        roi_type=meta["ROIType"],
    )
    roi = roi_crud.create(db, params=params)
    return roi


def _import_roi_point(db: Session, meta: dict, roi_id: int):
    params = ROIPointCreateModel(
        roi_id=roi_id,
        order_number=meta["OrderNumber"],
        slide_x_pos_um=meta["SlideXPosUm"],
        slide_y_pos_um=meta["SlideYPosUm"],
        panorama_pixel_x_pos=meta["PanoramaPixelXPos"],
        panorama_pixel_y_pos=meta["PanoramaPixelYPos"],
    )
    roi_point = roi_point_crud.create(db, params=params)
    return roi_point


def _import_acquisition(db: Session, meta: dict, roi_id: int):
    params = AcquisitionCreateModel(
        roi_id=roi_id,
        description=meta["Description"],
        order_number=meta["OrderNumber"],
        ablation_power=meta["AblationPower"],
        ablation_distance_between_shots_x=meta["AblationDistanceBetweenShotsX"],
        ablation_distance_between_shots_y=meta["AblationDistanceBetweenShotsY"],
        ablation_frequency=meta["AblationFrequency"],
        signal_type=meta["SignalType"],
        dual_count_start=meta["DualCountStart"],
        data_start_offset=meta["DataStartOffset"],
        data_end_offset=meta["DataEndOffset"],
        start_timestamp=meta["StartTimeStamp"],
        end_timestamp=meta["EndTimeStamp"],
        after_ablation_image_end_offset=meta["AfterAblationImageEndOffset"],
        after_ablation_image_start_offset=meta["AfterAblationImageStartOffset"],
        before_ablation_image_end_offset=meta["BeforeAblationImageEndOffset"],
        before_ablation_image_start_offset=meta["BeforeAblationImageStartOffset"],
        roi_start_x_pos_um=meta["ROIStartXPosUm"],
        roi_start_y_pos_um=meta["ROIStartYPosUm"],
        roi_end_x_pos_um=meta["ROIEndXPosUm"],
        roi_end_y_pos_um=meta["ROIEndYPosUm"],
        movement_type=meta["MovementType"],
        segment_data_format=meta["SegmentDataFormat"],
        value_bytes=meta["ValueBytes"],
        max_y=meta["MaxY"],
        max_x=meta["MaxX"],
        plume_start=meta["PlumeStart"],
        plume_end=meta["PlumeEnd"],
        template=meta["Template"],
        meta=meta,
    )
    acquisition = acquisition_crud.create(db, params=params)
    return acquisition


def _import_channel(db: Session, acquisition_id: int, index: int, imc_acquisition, acquisition_item):
    label = imc_acquisition.channel_labels[index]
    metal = imc_acquisition.channel_metals[index]
    mass = imc_acquisition.channel_mass[index]
    img = imc_acquisition.get_img_by_label(label)
    # metas = [a for a in acquisition_item.childs["AcquisitionChannel"].values() if label in a.properties["ChannelLabel"]]
    params = ChannelCreateModel(
        acquisition_id=acquisition_id,
        metal=metal,
        label=label,
        mass=mass,
        max_intensity=img.max(),
        min_intensity=img.min(),
        # meta=metas[0].properties if len(metas) > 0 else None
    )
    channel = channel_crud.create(db, params=params)
    np.save(os.path.join(channel.location, "origin.npy"), img)
    return channel


@timeit
def import_mcd(db: Session, uri: str, experiment_id: int):
    with McdParser(uri) as mcd:
        slide_item = mcd.meta.objects["Slide"]["0"]
        slide = _import_slide(db, slide_item.properties, experiment_id)

        for panorama_item in slide_item.childs["Panorama"].values():
            panorama = _import_panorama(db, panorama_item.properties, slide.id)

            if "AcquisitionROI" in panorama_item.childs:
                for roi_item in panorama_item.childs["AcquisitionROI"].values():
                    roi = _import_roi(db, roi_item.properties, panorama.id)

                    for roi_point_item in roi_item.childs["ROIPoint"].values():
                        roi_point = _import_roi_point(db, roi_point_item.properties, roi.id)

                    for acquisition_item in roi_item.childs["Acquisition"].values():
                        acquisition = _import_acquisition(db, acquisition_item.properties, roi.id)

                        imc_acquisition = mcd.get_imc_acquisition(acquisition_item.properties["ID"])
                        for i in range(imc_acquisition.n_channels):
                            channel = _import_channel(db, acquisition.id, i, imc_acquisition, acquisition_item)
