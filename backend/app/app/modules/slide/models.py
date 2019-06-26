from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel

from app.modules.panorama.models import PanoramaDatasetModel


# Properties to receive via API on creation
class SlideCreateModel(BaseModel):
    experiment_id: int
    uid: Optional[str]
    description: Optional[str]
    filename: Optional[str]
    slide_type: Optional[str]
    width_um: Optional[int]
    height_um: Optional[int]
    image_end_offset: Optional[int]
    image_start_offset: Optional[int]
    image_file: Optional[str]
    meta: Optional[dict]


# Shared properties
class SlideModel(BaseModel):
    id: int
    experiment_id: int
    uid: Optional[str]
    description: Optional[str]
    filename: Optional[str]
    slide_type: Optional[str]
    width_um: Optional[int]
    height_um: Optional[int]
    image_end_offset: Optional[int]
    image_start_offset: Optional[int]
    image_file: Optional[str]
    meta: Optional[dict]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


# Full slide dataset
class SlideDatasetModel(SlideModel):
    panoramas: List[PanoramaDatasetModel]

    class Config:
        orm_mode = True
