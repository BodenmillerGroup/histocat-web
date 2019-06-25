from datetime import datetime
from typing import Optional

from pydantic import BaseModel


# Properties to receive via API on creation
class PanoramaCreateModel(BaseModel):
    slide_id: int
    description: Optional[str]
    slide_x1_pos_um: Optional[float]
    slide_y1_pos_um: Optional[float]
    slide_x2_pos_um: Optional[float]
    slide_y2_pos_um: Optional[float]
    slide_x3_pos_um: Optional[float]
    slide_y3_pos_um: Optional[float]
    slide_x4_pos_um: Optional[float]
    slide_y4_pos_um: Optional[float]
    image_end_offset: Optional[int]
    image_start_offset: Optional[int]
    pixel_width: Optional[int]
    pixel_height: Optional[int]
    image_format: Optional[str]
    pixel_scale_coef: Optional[float]
    meta: Optional[dict]


# Shared properties
class PanoramaModel(BaseModel):
    id: int
    slide_id: int
    description: Optional[str]
    slide_x1_pos_um: Optional[float]
    slide_y1_pos_um: Optional[float]
    slide_x2_pos_um: Optional[float]
    slide_y2_pos_um: Optional[float]
    slide_x3_pos_um: Optional[float]
    slide_y3_pos_um: Optional[float]
    slide_x4_pos_um: Optional[float]
    slide_y4_pos_um: Optional[float]
    image_end_offset: Optional[int]
    image_start_offset: Optional[int]
    pixel_width: Optional[int]
    pixel_height: Optional[int]
    image_format: Optional[str]
    pixel_scale_coef: Optional[float]
    meta: Optional[dict]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True
