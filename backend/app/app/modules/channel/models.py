from datetime import datetime
from typing import List, Optional, Dict

from pydantic import BaseModel


# Properties to receive via API on creation
class ChannelCreateModel(BaseModel):
    acquisition_id: int
    origin_id: int
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
    origin_id: int
    metal: str
    label: str
    mass: int
    max_intensity: Optional[float]
    min_intensity: Optional[float]
    meta: Dict[str, Optional[str]]
    created_at: datetime

    class Config:
        orm_mode = True


class ChannelStatsModel(BaseModel):
    """
    Channel stats model
    """

    hist: List[int]
    edges: List[float]


class ChannelSettingsModel(BaseModel):
    id: int
    color: Optional[str]
    customLabel: Optional[str]
    min: Optional[float]
    max: Optional[float]


class FilterModel(BaseModel):
    apply: bool
    type: str
    settings: Optional[dict]


class LegendModel(BaseModel):
    apply: bool
    fontScale: float
    showIntensity: bool


class ScalebarModel(BaseModel):
    apply: bool
    settings: Optional[dict]


class MaskSettingsModel(BaseModel):
    apply: bool
    colorize: Optional[bool]
    location: str
    settings: Optional[dict]


class ChannelStackModel(BaseModel):
    datasetId: Optional[int]
    filter: FilterModel
    legend: LegendModel
    scalebar: ScalebarModel
    channels: List[ChannelSettingsModel]
    mask: Optional[MaskSettingsModel]
    format: Optional[str] = "png"
