from datetime import datetime
from typing import Optional, List, Dict

from pydantic import BaseModel

from app.modules.roi.models import ROIDatasetModel


# Properties to receive via API on creation
class PanoramaCreateModel(BaseModel):
    slide_id: int
    metaname: str
    original_id: int
    meta: Dict[str, Optional[str]]


# Shared properties
class PanoramaModel(BaseModel):
    id: int
    slide_id: int
    metaname: str
    original_id: int
    meta: Dict[str, Optional[str]]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


# Full panorama dataset
class PanoramaDatasetModel(PanoramaModel):
    rois: List[ROIDatasetModel]

    class Config:
        orm_mode = True
