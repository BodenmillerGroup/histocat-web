from datetime import datetime
from typing import Any, Optional, Dict, List

from histocat.modules.analysis.dto import PlotSeriesDto
from pydantic import BaseModel


class ResultCreateDto(BaseModel):
    dataset_id: int
    status: str
    name: Optional[str]
    description: Optional[str]
    pipeline: Any
    input: Any


class ResultUpdateDto(BaseModel):
    name: Optional[str]
    description: Optional[str]
    output: Optional[Any]


class ResultDto(BaseModel):
    id: int
    dataset_id: int
    status: str
    name: Optional[str]
    description: Optional[str]
    pipeline: Any
    input: Any
    output: Any
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


class MappingDto(BaseModel):
    """2D data mapping"""
    x: PlotSeriesDto
    y: PlotSeriesDto


class ColorsDto(BaseModel):
    """Color data mapping"""
    type: str
    name: str
    data: List[Any]


class ResultDataDto(BaseModel):
    cellIds: List[str]
    acquisitionIds: List[int]
    objectNumbers: List[int]
    markers: List[str]
    x: List[float]
    y: List[float]
    mappings: Optional[Dict[str, MappingDto]]


class ColorsDataDto(BaseModel):
    cellIds: List[str]
    colors: ColorsDto
