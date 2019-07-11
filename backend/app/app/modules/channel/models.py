from datetime import datetime
from typing import List, Optional, Dict

from pydantic import BaseModel


# Properties to receive via API on creation
class ChannelCreateModel(BaseModel):
    acquisition_id: int
    metaname: str
    original_id: int
    metal: str
    label: str
    mass: int
    max_intensity: Optional[int]
    min_intensity: Optional[int]
    meta: Dict[str, str]


# Shared properties
class ChannelModel(BaseModel):
    id: int
    acquisition_id: int
    metaname: str
    original_id: int
    metal: str
    label: str
    mass: int
    max_intensity: Optional[int]
    min_intensity: Optional[int]
    meta: Dict[str, str]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


# Channel stats model
class ChannelStatsModel(BaseModel):
    hist: List[int]
    bins: List[float]
