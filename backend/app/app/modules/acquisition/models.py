from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel

from app.modules.channel.models import ChannelModel


# Properties to receive via API on creation
class AcquisitionCreateModel(BaseModel):
    roi_id: int
    description: Optional[str]
    order_number: Optional[int]
    ablation_power: Optional[float]
    ablation_distance_between_shots_x: Optional[float]
    ablation_distance_between_shots_y: Optional[float]
    ablation_frequency: Optional[float]
    signal_type: Optional[str]
    dual_count_start: Optional[int]
    data_start_offset: Optional[int]
    data_end_offset: Optional[int]
    start_timestamp: Optional[datetime]
    end_timestamp: Optional[datetime]
    after_ablation_image_end_offset: Optional[int]
    after_ablation_image_start_offset: Optional[int]
    before_ablation_image_end_offset: Optional[int]
    before_ablation_image_start_offset: Optional[int]
    roi_start_x_pos_um: Optional[float]
    roi_start_y_pos_um: Optional[float]
    roi_end_x_pos_um: Optional[float]
    roi_end_y_pos_um: Optional[float]
    movement_type: Optional[str]
    segment_data_format: Optional[str]
    value_bytes: Optional[int]
    max_y: Optional[int]
    max_x: Optional[int]
    plume_start: Optional[int]
    plume_end: Optional[int]
    template: Optional[str]
    meta: Optional[dict]


# Shared properties
class AcquisitionModel(BaseModel):
    id: int
    roi_id: int
    description: Optional[str]
    order_number: Optional[int]
    ablation_power: Optional[float]
    ablation_distance_between_shots_x: Optional[float]
    ablation_distance_between_shots_y: Optional[float]
    ablation_frequency: Optional[float]
    signal_type: Optional[str]
    dual_count_start: Optional[int]
    data_start_offset: Optional[int]
    data_end_offset: Optional[int]
    start_timestamp: Optional[datetime]
    end_timestamp: Optional[datetime]
    after_ablation_image_end_offset: Optional[int]
    after_ablation_image_start_offset: Optional[int]
    before_ablation_image_end_offset: Optional[int]
    before_ablation_image_start_offset: Optional[int]
    roi_start_x_pos_um: Optional[float]
    roi_start_y_pos_um: Optional[float]
    roi_end_x_pos_um: Optional[float]
    roi_end_y_pos_um: Optional[float]
    movement_type: Optional[str]
    segment_data_format: Optional[str]
    value_bytes: Optional[int]
    max_y: Optional[int]
    max_x: Optional[int]
    plume_start: Optional[int]
    plume_end: Optional[int]
    template: Optional[str]
    meta: Optional[dict]
    created_at: datetime

    class Config:
        orm_mode = True


# Acquisition dataset with children items
class AcquisitionDatasetModel(AcquisitionModel):
    channels: Optional[List[ChannelModel]]

    class Config:
        orm_mode = True
