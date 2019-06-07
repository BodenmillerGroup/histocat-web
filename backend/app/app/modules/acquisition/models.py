from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel

from app.modules.channel.models import ChannelModel


# Properties to receive via API on creation
class AcquisitionCreateModel(BaseModel):
    slide_id: int
    name: str
    width: int
    height: int
    description: Optional[str] = None
    meta: Optional[dict] = None


# Shared properties
class AcquisitionModel(BaseModel):
    id: int
    slide_id: int
    name: str
    width: int
    height: int
    location: Optional[str]
    description: Optional[str]
    meta: Optional[dict]
    created_at: datetime


# Acquisition dataset with children items
class AcquisitionDatasetModel(AcquisitionModel):
    channels: Optional[List[ChannelModel]]
