from datetime import datetime
from typing import Optional

from pydantic import BaseModel


# Properties to receive via API on creation
class ROIPointCreateModel(BaseModel):
    roi_id: int
    order_number: Optional[int]
    slide_x_pos_um: Optional[float]
    slide_y_pos_um: Optional[float]
    panorama_pixel_x_pos: Optional[int]
    panorama_pixel_y_pos: Optional[int]


# Shared properties
class ROIPointModel(BaseModel):
    id: int
    roi_id: int
    order_number: Optional[int]
    slide_x_pos_um: Optional[float]
    slide_y_pos_um: Optional[float]
    panorama_pixel_x_pos: Optional[int]
    panorama_pixel_y_pos: Optional[int]
    created_at: datetime

    class Config:
        orm_mode = True
