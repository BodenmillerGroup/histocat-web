from datetime import datetime
from typing import List, Optional, Dict

from pydantic import BaseModel

from app.modules.channel.models import ChannelModel


# Properties to receive via API on creation
class AcquisitionCreateModel(BaseModel):
    roi_id: int
    metaname: str
    original_id: int
    meta: Dict[str, Optional[str]]


# Shared properties
class AcquisitionModel(BaseModel):
    id: int
    roi_id: int
    metaname: str
    original_id: int
    meta: Dict[str, Optional[str]]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


# Acquisition dataset with children items
class AcquisitionDatasetModel(AcquisitionModel):
    channels: Optional[List[ChannelModel]]

    class Config:
        orm_mode = True
