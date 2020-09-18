from datetime import datetime
from typing import Optional, Sequence

from pydantic import BaseModel

from histocat.modules.slide.dto import SlideDatasetDto


class ProjectCreateDto(BaseModel):
    name: str
    description: Optional[str]
    tags: Optional[Sequence[str]]


class ProjectUpdateDto(BaseModel):
    name: str
    description: Optional[str]
    tags: Optional[Sequence[str]]


class ProjectDto(BaseModel):
    id: int
    group_id: int
    name: Optional[str]
    description: Optional[str]
    tags: Optional[Sequence[str]]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


class ProjectFullDto(ProjectDto):
    """Full set of project data."""

    slides: Sequence[SlideDatasetDto]

    class Config:
        orm_mode = True
