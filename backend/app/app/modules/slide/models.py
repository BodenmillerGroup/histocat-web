from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel

from app.modules.acquisition.models import AcquisitionDatasetModel


# Properties to receive via API on creation
class SlideCreateModel(BaseModel):
    experiment_id: int
    name: str
    filename: Optional[str] = None
    width_um: Optional[int] = None
    height_um: Optional[int] = None
    description: Optional[str] = None
    meta: Optional[dict] = None


# Shared properties
class SlideModel(BaseModel):
    id: int
    experiment_id: int
    name: str
    filename: Optional[str]
    width_um: Optional[int]
    height_um: Optional[int]
    location: Optional[str]
    description: Optional[str]
    meta: Optional[dict]
    created_at: datetime

    class Config:
        orm_mode = True


# Full slide dataset
class SlideDatasetModel(SlideModel):
    acquisitions: Optional[List[AcquisitionDatasetModel]]

    class Config:
        orm_mode = True
