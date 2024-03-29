from datetime import datetime
from typing import Any, Dict, List, Optional

from pydantic import BaseModel

from histocat.core.analysis.dto import PlotSeriesDto


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
    status: Optional[str]


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
    data: List[str]


class ResultDataDto(BaseModel):
    cellIds: List[str]
    markers: List[str]
    mappings: Optional[Dict[str, MappingDto]]


class ColorsDataDto(BaseModel):
    cellIds: List[str]
    colors: ColorsDto
