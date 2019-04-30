from datetime import datetime
from typing import Optional, List

from pydantic import BaseModel

from app.modules.acquisition.models import AcquisitionDatasetModel


# Shared properties
class SlideBaseModel(BaseModel):
    id: int = None
    experiment_id: int = None
    name: Optional[str] = None
    filename: Optional[str] = None
    width_um: Optional[int] = None
    height_um: Optional[int] = None
    location: Optional[str] = None
    description: Optional[str] = None
    meta: Optional[dict] = None
    created_at: datetime


# Properties to receive via API on creation
class SlideCreateModel(BaseModel):
    experiment_id: int
    name: str
    filename: Optional[str] = None
    width_um: Optional[int] = None
    height_um: Optional[int] = None
    description: Optional[str] = None
    meta: Optional[dict] = None


# Properties to receive via API on update
class SlideUpdateModel(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    meta: Optional[dict] = None


# Additional properties to return via API
class SlideModel(SlideBaseModel):
    pass


# Additional properties stored in DB
class SlideInDBModel(SlideBaseModel):
    pass


# Full slide dataset
class SlideDatasetModel(SlideModel):
    acquisitions: Optional[List[AcquisitionDatasetModel]]
