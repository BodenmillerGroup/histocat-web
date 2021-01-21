from datetime import datetime
from typing import Optional, Sequence

from pydantic import BaseModel

from histocat.core.member.dto import MemberDto


class GroupCreateDto(BaseModel):
    name: str
    description: Optional[str]
    url: Optional[str]
    is_open: bool
    tags: Optional[Sequence[str]]


class GroupUpdateDto(BaseModel):
    name: str
    description: Optional[str]
    url: Optional[str]
    is_open: bool
    tags: Optional[Sequence[str]]


class GroupDto(BaseModel):
    id: int
    name: str
    description: Optional[str]
    url: Optional[str]
    is_open: bool
    tags: Optional[Sequence[str]]
    location: Optional[str]
    created_at: datetime

    members: Sequence[MemberDto]

    class Config:
        orm_mode = True
