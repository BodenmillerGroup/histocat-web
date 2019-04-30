from datetime import datetime
from typing import Optional, List

from pydantic import BaseModel

from app.modules.channel.models import ChannelModel


# Shared properties
class AcquisitionBaseModel(BaseModel):
    id: int = None
    slide_id: int = None
    name: Optional[str] = None
    width: Optional[int] = None
    height: Optional[int] = None
    location: Optional[str] = None
    description: Optional[str] = None
    meta: Optional[dict] = None
    created_at: datetime


# Properties to receive via API on creation
class AcquisitionCreateModel(BaseModel):
    slide_id: int
    name: str
    width: int
    height: int
    description: Optional[str] = None
    meta: Optional[dict] = None


# Properties to receive via API on update
class AcquisitionUpdateModel(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    meta: Optional[dict] = None


# Additional properties to return via API
class AcquisitionModel(AcquisitionBaseModel):
    pass


# Additional properties stored in DB
class AcquisitionInDBModel(AcquisitionBaseModel):
    pass


# Acquisition dataset with children items
class AcquisitionDatasetModel(AcquisitionModel):
    channels: Optional[List[ChannelModel]]
