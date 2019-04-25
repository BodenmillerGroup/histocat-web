from typing import Optional

from datetime import datetime
from pydantic import BaseModel


# Shared properties
class SlideBaseModel(BaseModel):
    id: int = None
    experiment_id: int = None
    name: Optional[str] = None
    filename: Optional[str] = None
    width_um: Optional[int] = None
    height_um: Optional[int] = None
    description: Optional[str] = None
    meta: Optional[dict] = None


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
    location: Optional[str] = None


# Additional properties stored in DB
class SlideInDBModel(SlideBaseModel):
    created_at: datetime
