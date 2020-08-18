from datetime import datetime
from typing import Optional, Sequence

from pydantic import BaseModel

from histocat.modules.slide.dto import SlideDatasetDto


class ExperimentCreateDto(BaseModel):
    name: str
    description: Optional[str]
    meta: Optional[dict]
    tags: Optional[Sequence[str]]
    is_public: Optional[bool]


class ExperimentUpdateDto(BaseModel):
    name: str
    description: Optional[str]
    meta: Optional[dict]
    tags: Optional[Sequence[str]]
    is_public: Optional[bool]


class ExperimentDto(BaseModel):
    id: int
    group_id: int
    member_id: int
    name: Optional[str]
    description: Optional[str]
    meta: Optional[dict]
    tags: Optional[Sequence[str]]
    is_public: bool
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


class ExperimentDatasetDto(ExperimentDto):
    """Full experiment dataset."""

    slides: Sequence[SlideDatasetDto]

    class Config:
        orm_mode = True
