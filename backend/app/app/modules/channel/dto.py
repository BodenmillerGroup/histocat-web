from typing import Dict, List, Optional

from pydantic import BaseModel


class ChannelCreateDto(BaseModel):
    acquisition_id: int
    origin_id: int
    metal: str
    label: str
    mass: int
    max_intensity: Optional[float]
    min_intensity: Optional[float]
    meta: Dict[str, Optional[str]]


class ChannelDto(BaseModel):
    id: int
    acquisition_id: int
    origin_id: int
    metal: str
    label: str
    mass: int
    max_intensity: Optional[float]
    min_intensity: Optional[float]
    meta: Dict[str, Optional[str]]

    class Config:
        orm_mode = True


class ChannelStatsDto(BaseModel):
    """Channel stats model."""
    hist: List[int]
    edges: List[float]


class ChannelSettingsDto(BaseModel):
    id: int
    color: Optional[str]
    customLabel: Optional[str]
    min: Optional[float]
    max: Optional[float]


class FilterDto(BaseModel):
    apply: bool
    type: str
    settings: Optional[dict]


class LegendDto(BaseModel):
    apply: bool
    fontScale: float
    showIntensity: bool


class ScalebarDto(BaseModel):
    apply: bool
    settings: Optional[dict]


class MaskSettingsDto(BaseModel):
    apply: bool
    colorize: Optional[bool]
    location: str
    settings: Optional[dict]


class ChannelStackDto(BaseModel):
    datasetId: Optional[int]
    filter: FilterDto
    legend: LegendDto
    scalebar: ScalebarDto
    channels: List[ChannelSettingsDto]
    mask: Optional[MaskSettingsDto]
    format: Optional[str] = "png"
