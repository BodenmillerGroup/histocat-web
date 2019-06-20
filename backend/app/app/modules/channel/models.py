from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel


# Properties to receive via API on creation
class ChannelCreateModel(BaseModel):
    acquisition_id: int
    name: str
    metal: str
    mass: int
    max_intensity: int
    min_intensity: int
    meta: Optional[dict] = None


# Shared properties
class ChannelModel(BaseModel):
    id: int
    acquisition_id: int
    name: str
    metal: str
    mass: int
    max_intensity: int
    min_intensity: int
    location: Optional[str]
    meta: Optional[dict]
    created_at: datetime

    class Config:
        orm_mode = True


# Channel stats model
class ChannelStatsModel(BaseModel):
    hist: List[int]
    bins: List[float]
