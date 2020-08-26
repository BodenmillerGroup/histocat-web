from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Dict

from imctools import data
from sqlalchemy.orm import Session

from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.io.utils import copy_dir
from histocat.modules.acquisition import service as acquisition_service
from histocat.modules.acquisition.dto import AcquisitionCreateDto
from histocat.modules.acquisition.models import AcquisitionModel
from histocat.modules.panorama import service as panorama_service
from histocat.modules.panorama.dto import PanoramaCreateDto
from histocat.modules.panorama.models import PanoramaModel
from histocat.modules.slide import service as slide_service
from histocat.modules.slide.dto import SlideCreateDto
from histocat.modules.slide.models import SlideModel

logger = logging.getLogger(__name__)


def import_imcfolder(db: Session, session_filename: str, experiment_id: int):
    """
    Import slides from the folder compatible with IMC pipeline
    """

    session_path = Path(session_filename)
    src_folder = session_path.parent
    session = data.Session.load(session_filename)
    basename = session.name  # schema_path.stem.replace("_schema", "")

    slide_map: Dict[int, SlideModel] = dict()
    for slide_key, slide_item in session.slides.items():
        if slide_key is None:
            continue
        item = slide_service.get_by_name(db, experiment_id=experiment_id, name=basename)
        if item:
            logger.warning(f"The slide with the name [{basename}] already exists in the experiment [{experiment_id}]")
            return

        slide = _import_slide(db, slide_item, experiment_id, basename)
        slide_map[slide.origin_id] = slide

        origin_location = os.path.join(slide.location, "origin")
        copy_dir(src_folder, origin_location)

    panorama_map: Dict[int, PanoramaModel] = dict()
    for panorama_key, panorama_item in session.panoramas.items():
        if panorama_key is None:
            continue
        slide = slide_map.get(panorama_item.slide_id)
        panorama = _import_panorama(db, panorama_item, slide)
        panorama_map[panorama.origin_id] = panorama

    acquisition_map: Dict[int, AcquisitionModel] = dict()
    for acquisition_key, acquisition_item in session.acquisitions.items():
        if acquisition_key is None:
            continue
        slide = slide_map.get(acquisition_item.slide_id)
        acquisition = _import_acquisition(db, acquisition_item, slide, basename)
        acquisition_map[acquisition.origin_id] = acquisition

    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(experiment_id, "slide_imported"))


def _import_slide(db: Session, item: data.Slide, experiment_id: int, name: str):
    params = SlideCreateDto(
        experiment_id=experiment_id,
        name=name,
        origin_id=item.id,
        width_um=item.width_um,
        height_um=item.height_um,
        has_slide_image=item.has_slide_image,
        meta=item.metadata,
        session_meta=item.session.metadata,
    )
    slide = slide_service.create(db, params=params)
    return slide


def _import_panorama(db: Session, item: data.Panorama, slide: SlideModel):
    origin_location = os.path.join(slide.location, "origin")
    location = os.path.join(origin_location, f"{item.metaname}_pano.png",)
    params = PanoramaCreateDto(
        slide_id=slide.id,
        origin_id=item.id,
        image_type=item.image_type,
        description=item.description,
        x1=item.x1,
        y1=item.y1,
        x2=item.x2,
        y2=item.y2,
        x3=item.x3,
        y3=item.y3,
        x4=item.x4,
        y4=item.y4,
        rotation_angle=item.rotation_angle,
        location=location if item.image_type != "Default" else None,
        meta=item.metadata,
    )
    panorama = panorama_service.create(db, params=params)
    return panorama


def _import_acquisition(db: Session, item: data.Acquisition, slide: SlideModel, basename: str):
    origin_location = os.path.join(slide.location, "origin")
    location = os.path.join(origin_location, f"{item.metaname}_ac.ome.tiff",)

    # Import channel data
    channels = dict()
    for channel_key, channel_item in item.channels.items():
        channel_data = channel_item.__getstate__()
        del channel_data["metadata"]
        channels[channel_item.name] = channel_data

    params = AcquisitionCreateDto(
        slide_id=slide.id,
        origin_id=item.id,
        description=item.description,
        max_x=item.max_x,
        max_y=item.max_y,
        signal_type=item.signal_type,
        segment_data_format=item.segment_data_format,
        ablation_frequency=item.ablation_frequency,
        ablation_power=item.ablation_power,
        start_timestamp=item.start_timestamp.isoformat() if item.is_valid else None,
        end_timestamp=item.end_timestamp.isoformat() if item.is_valid else None,
        movement_type=item.movement_type,
        ablation_distance_between_shots_x=item.ablation_distance_between_shots_x,
        ablation_distance_between_shots_y=item.ablation_distance_between_shots_y,
        template=item.template,
        roi_start_x_pos_um=item.roi_start_x_pos_um,
        roi_start_y_pos_um=item.roi_start_y_pos_um,
        roi_end_x_pos_um=item.roi_end_x_pos_um,
        roi_end_y_pos_um=item.roi_end_y_pos_um,
        has_before_ablation_image=item.has_before_ablation_image,
        has_after_ablation_image=item.has_after_ablation_image,
        is_valid=item.is_valid,
        channels=channels,
        location=location if item.is_valid else None,
        meta=item.metadata,
    )
    acquisition = acquisition_service.create(db, params)
    return acquisition
