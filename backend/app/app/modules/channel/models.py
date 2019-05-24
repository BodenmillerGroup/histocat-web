from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel


# Shared properties
class ChannelBaseModel(BaseModel):
    id: int = None
    acquisition_id: int = None
    name: Optional[str] = None
    metal: Optional[str] = None
    mass: Optional[int] = None
    max_intensity: Optional[int] = None
    min_intensity: Optional[int] = None
    location: Optional[str] = None
    meta: Optional[dict] = None
    created_at: datetime


# Properties to receive via API on creation
class ChannelCreateModel(BaseModel):
    acquisition_id: int
    name: str
    metal: str
    mass: int
    max_intensity: int
    min_intensity: int
    meta: Optional[dict] = None


# Properties to receive via API on update
class ChannelUpdateModel(BaseModel):
    name: Optional[str] = None
    meta: Optional[dict] = None


# Additional properties to return via API
class ChannelModel(ChannelBaseModel):
    pass


# Additional properties stored in DB
class ChannelInDBModel(ChannelBaseModel):
    pass


# Channel stats model
class ChannelStatsModel(BaseModel):
    hist: List[int]
    bins: List[float]
