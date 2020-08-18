from datetime import datetime
from typing import Optional

from pydantic import BaseModel


class GroupCreateDto(BaseModel):
    name: str
    description: Optional[str]
    url: Optional[str]
    is_open: bool


class GroupUpdateDto(BaseModel):
    name: str
    description: Optional[str]
    url: Optional[str]
    is_open: bool


class GroupDto(BaseModel):
    id: int
    name: str
    description: Optional[str]
    url: Optional[str]
    is_open: bool
    created_at: datetime

    class Config:
        orm_mode = True
