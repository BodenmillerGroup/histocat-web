import logging
import os
from pathlib import Path
from typing import Dict, List

from imctools.io.ometiff.ometiffparser import OmeTiffParser
from sqlalchemy.orm import Session

from histocat.core.acquisition import service as acquisition_service
from histocat.core.acquisition.dto import AcquisitionCreateDto
from histocat.core.acquisition.models import AcquisitionModel
from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.core.slide import service as slide_service
from histocat.core.slide.dto import SlideCreateDto
from histocat.core.slide.models import SlideModel
from histocat.worker.io.utils import copy_dir

logger = logging.getLogger(__name__)


def import_ometifffolder(db: Session, slide_name: str, ome_tiff_files: List[str], project_id: int):
    """
    Import slide from the folder with OME-TIFFs
    """

    item = slide_service.get_by_name(db, project_id=project_id, name=slide_name)
    if item:
        logger.warning(f"The slide with the name [{slide_name}] already exists in the project [{project_id}]")
        return

    src_folder = Path(ome_tiff_files[0]).parent
    slide = _import_slide(db, project_id, slide_name)
    origin_location = os.path.join(slide.location, "origin")
    copy_dir(src_folder, origin_location)

    acquisition_map: Dict[int, AcquisitionModel] = dict()
    for ome_tiff_file in ome_tiff_files:
        acquisition = _import_acquisition(db, slide, ome_tiff_file)
        acquisition_map[acquisition.origin_id] = acquisition

    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(project_id, "slide_imported"))


def _import_slide(db: Session, project_id: int, slide_name: str):
    params = SlideCreateDto(
        project_id=project_id,
        name=slide_name,
        origin_id=0,
        width_um=0,
        height_um=0,
        has_slide_image=False,
        meta={},
        session_meta={},
    )
    slide = slide_service.create(db, params=params)
    return slide


def _import_acquisition(db: Session, slide: SlideModel, ometiff_filename: str):
    origin_location = os.path.join(slide.location, "origin")
    location = os.path.join(origin_location, Path(ometiff_filename).name)
    ometiff = OmeTiffParser(location, slide.id)
    ac_data = ometiff.get_acquisition_data()
    ac = ac_data.acquisition

    # Import channel data
    channels = dict()
    for key, channel in ac.channels.items():
        channel_data = channel.__getstate__()
        del channel_data["metadata"]
        channel_data["customLabel"] = channel.label
        channel_data["mass"] = int("".join([n for n in channel.name if n.isdigit()]))
        channels[channel.name] = channel_data

    params = AcquisitionCreateDto(
        slide_id=slide.id,
        origin_id=ac.id,
        description=ac.description,
        max_x=ac.max_x,
        max_y=ac.max_y,
        signal_type=ac.signal_type,
        segment_data_format=ac.segment_data_format,
        ablation_frequency=ac.ablation_frequency,
        ablation_power=ac.ablation_power,
        start_timestamp=ac.start_timestamp.isoformat() if ac.start_timestamp else None,
        end_timestamp=ac.end_timestamp.isoformat() if ac.start_timestamp else None,
        movement_type=ac.movement_type,
        ablation_distance_between_shots_x=ac.ablation_distance_between_shots_x,
        ablation_distance_between_shots_y=ac.ablation_distance_between_shots_y,
        template=ac.template,
        roi_start_x_pos_um=ac.roi_start_x_pos_um,
        roi_start_y_pos_um=ac.roi_start_y_pos_um,
        roi_end_x_pos_um=ac.roi_end_x_pos_um,
        roi_end_y_pos_um=ac.roi_end_y_pos_um,
        has_before_ablation_image=ac.has_before_ablation_image,
        has_after_ablation_image=ac.has_after_ablation_image,
        is_valid=ac_data.is_valid,
        channels=channels,
        location=location if ac_data.is_valid else None,
        meta=ac.metadata,
    )
    acquisition = acquisition_service.create(db, params)
    return acquisition
