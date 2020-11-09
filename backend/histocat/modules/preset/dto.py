from datetime import datetime
from typing import Any, Optional

from pydantic import BaseModel


class PresetCreateDto(BaseModel):
    project_id: int
    name: Optional[str]
    description: Optional[str]
    data: Optional[Any]


class PresetDto(BaseModel):
    id: int
    project_id: Optional[int]
    name: Optional[str]
    description: Optional[str]
    data: Optional[Any]
    created_at: datetime

    class Config:
        orm_mode = True
