from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel


# Properties to receive via API on creation
class ChannelCreateModel(BaseModel):
    acquisition_id: int
    channel_name: Optional[str]
    channel_label: Optional[str]
    order_number: Optional[int]
    mass: Optional[int]
    max_intensity: Optional[int]
    min_intensity: Optional[int]
    meta: Optional[dict]


# Shared properties
class ChannelModel(BaseModel):
    id: int
    acquisition_id: int
    channel_name: Optional[str]
    channel_label: Optional[str]
    order_number: Optional[int]
    mass: Optional[int]
    max_intensity: Optional[int]
    min_intensity: Optional[int]
    meta: Optional[dict]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


# Channel stats model
class ChannelStatsModel(BaseModel):
    hist: List[int]
    bins: List[float]
