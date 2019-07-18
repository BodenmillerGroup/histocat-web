from datetime import datetime
from typing import List, Optional, Dict, Any

from pydantic import BaseModel


# Properties to receive via API on creation
class ChannelCreateModel(BaseModel):
    acquisition_id: int
    metaname: str
    original_id: int
    metal: str
    label: str
    mass: int
    max_intensity: Optional[float]
    min_intensity: Optional[float]
    meta: Dict[str, Optional[str]]


# Shared properties
class ChannelModel(BaseModel):
    id: int
    acquisition_id: int
    metaname: str
    original_id: int
    metal: str
    label: str
    mass: int
    max_intensity: Optional[float]
    min_intensity: Optional[float]
    meta: Dict[str, Optional[str]]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


# Channel stats model
class ChannelStatsModel(BaseModel):
    hist: List[int]
    edges: List[float]


class ChannelSettingsModel(BaseModel):
    id: int
    color: Optional[str]
    min: Optional[float]
    max: Optional[float]


class FilterModel(BaseModel):
    apply: bool
    type: str
    settings: Optional[dict]


class ChannelStackModel(BaseModel):
    filter: FilterModel
    channels: List[ChannelSettingsModel]
    format: Optional[str] = 'png'
