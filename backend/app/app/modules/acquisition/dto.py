from datetime import datetime
from typing import Dict, Optional, Any

from pydantic import BaseModel


class AcquisitionCreateDto(BaseModel):
    slide_id: int
    origin_id: int
    description: Optional[str]
    max_x: Optional[int]
    max_y: Optional[int]
    signal_type: Optional[str]
    segment_data_format: Optional[str]
    ablation_frequency: Optional[float]
    ablation_power: Optional[float]
    start_timestamp: Optional[str]
    end_timestamp: Optional[str]
    movement_type: Optional[str]
    ablation_distance_between_shots_x: Optional[float]
    ablation_distance_between_shots_y: Optional[float]
    template: Optional[str]
    roi_start_x_pos_um: Optional[float]
    roi_start_y_pos_um: Optional[float]
    roi_end_x_pos_um: Optional[float]
    roi_end_y_pos_um: Optional[float]
    has_before_ablation_image: Optional[bool]
    has_after_ablation_image: Optional[bool]
    is_valid: Optional[bool]
    location: Optional[str]
    meta: Dict[str, Optional[str]]
    channels: Optional[Dict[Any, Any]]


class AcquisitionDto(BaseModel):
    id: int
    slide_id: int
    origin_id: int
    description: Optional[str]
    max_x: Optional[int]
    max_y: Optional[int]
    signal_type: Optional[str]
    segment_data_format: Optional[str]
    ablation_frequency: Optional[float]
    ablation_power: Optional[float]
    start_timestamp: Optional[datetime]
    end_timestamp: Optional[datetime]
    movement_type: Optional[str]
    ablation_distance_between_shots_x: Optional[float]
    ablation_distance_between_shots_y: Optional[float]
    template: Optional[str]
    roi_start_x_pos_um: Optional[float]
    roi_start_y_pos_um: Optional[float]
    roi_end_x_pos_um: Optional[float]
    roi_end_y_pos_um: Optional[float]
    has_before_ablation_image: Optional[bool]
    has_after_ablation_image: Optional[bool]
    is_valid: Optional[bool]
    meta: Dict[str, Optional[str]]
    channels: Optional[Dict[Any, Any]]

    class Config:
        orm_mode = True


class AcquisitionDatasetDto(AcquisitionDto):
    class Config:
        orm_mode = True
