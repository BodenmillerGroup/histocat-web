from typing import Optional

from datetime import datetime
from pydantic import BaseModel


# Shared properties
class SlideBaseModel(BaseModel):
    id: int = None
    experiment_id: int = None
    name: Optional[str] = None
    description: Optional[str] = None
    location: Optional[str] = None
    meta: Optional[dict] = None


# Properties to receive via API on creation
class SlideInCreateModel(BaseModel):
    name: str
    experiment_id: int
    description: Optional[str] = None
    meta: Optional[dict] = None


# Properties to receive via API on update
class SlideInUpdateModel(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    meta: Optional[dict] = None


# Additional properties to return via API
class SlideModel(SlideBaseModel):
    pass


# Additional properties stored in DB
class SlideInDBModel(SlideBaseModel):
    created_at: datetime
