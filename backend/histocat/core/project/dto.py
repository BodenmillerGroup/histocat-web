from datetime import datetime
from typing import List, Optional, Sequence

from pydantic import BaseModel

from histocat.core.member.dto import MemberDto
from histocat.core.slide.dto import SlideDatasetDto


class ProjectCreateDto(BaseModel):
    name: str
    description: Optional[str]
    tags: Optional[List[str]]


class ProjectUpdateDto(BaseModel):
    name: str
    description: Optional[str]
    tags: Optional[List[str]]


class ProjectDto(BaseModel):
    id: int
    group_id: int
    name: Optional[str]
    description: Optional[str]
    tags: Optional[List[str]]
    location: Optional[str]
    created_at: datetime

    member: MemberDto

    class Config:
        orm_mode = True


class ProjectFullDto(ProjectDto):
    """Full set of project data."""

    slides: Sequence[SlideDatasetDto]

    class Config:
        orm_mode = True
