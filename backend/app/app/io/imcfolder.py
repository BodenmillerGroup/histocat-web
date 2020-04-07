from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Dict

from sqlalchemy.orm import Session

from imctools import data

from app.core.notifier import Message
from app.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from app.io.utils import copy_dir
from app.modules.acquisition import service as acquisition_service
from app.modules.acquisition.models import Acquisition
from app.modules.acquisition.dto import AcquisitionCreateDto
from app.modules.panorama import service as panorama_service
from app.modules.panorama.models import Panorama
from app.modules.panorama.dto import PanoramaCreateDto
from app.modules.slide import service as slide_service
from app.modules.slide.models import Slide
from app.modules.slide.dto import SlideCreateDto

logger = logging.getLogger(__name__)


def import_imcfolder(db: Session, session_filename: str, experiment_id: int, user_id: int):
    """
    Import slides from the folder compatible with IMC pipeline
    """

    session_path = Path(session_filename)
    src_folder = session_path.parent
    session = data.Session.load(session_filename)
    basename = session.name # schema_path.stem.replace("_schema", "")

    slide_map: Dict[int, Slide] = dict()
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

    panorama_map: Dict[int, Panorama] = dict()
    for panorama_key, panorama_item in session.panoramas.items():
        if panorama_key is None:
            continue
        slide = slide_map.get(panorama_item.slide_id)
        panorama = _import_panorama(db, panorama_item, slide)
        panorama_map[panorama.origin_id] = panorama

    acquisition_map: Dict[int, Acquisition] = dict()
    for acquisition_key, acquisition_item in session.acquisitions.items():
        if acquisition_key is None:
            continue
        slide = slide_map.get(acquisition_item.slide_id)
        acquisition = _import_acquisition(db, acquisition_item, slide, basename)
        acquisition_map[acquisition.origin_id] = acquisition

    # for channel_key, channel_item in channel_data.items():
    #     if channel_key is None:
    #         continue
    #     acquisition = acquisition_map.get(channel_item.get(mcdxmlparser.ACQUISITIONID))
    #     _import_channel(db, channel_item, acquisition)

    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(experiment_id, "slide_imported"))


def _import_slide(db: Session, slide_item: data.Slide, experiment_id: int, name: str):
    params = SlideCreateDto(
        experiment_id=experiment_id,
        name=name,
        origin_id=slide_item.id,
        meta=slide_item.metadata,
        xml_meta=None,
    )
    slide = slide_service.create(db, params=params)
    return slide


def _import_panorama(db: Session, panorama_item: data.Panorama, slide: Slide):
    params = PanoramaCreateDto(
        slide_id=slide.id,
        origin_id=panorama_item.id,
        meta=panorama_item.metadata
    )
    panorama = panorama_service.create(db, params=params)
    return panorama


def _import_acquisition(db: Session, acquisition_item: data.Acquisition, slide: Slide, basename: str):
    origin_location = os.path.join(slide.location, "origin")
    location = os.path.join(
        origin_location,
        f"{acquisition_item.meta_name}_ac.ome.tiff",
    )
    params = AcquisitionCreateDto(
        slide_id=slide.id,
        origin_id=acquisition_item.id,
        location=location,
        meta=acquisition_item.metadata
    )
    acquisition = acquisition_service.create(db, params)
    return acquisition
