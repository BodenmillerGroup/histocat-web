from datetime import datetime
from typing import Any, Dict, Optional

from pydantic import BaseModel


class PresetCreateDto(BaseModel):
    experiment_id: int
    name: Optional[str]
    description: Optional[str]
    data: Optional[Any]


class PresetDto(BaseModel):
    id: int
    experiment_id: int
    name: Optional[str]
    description: Optional[str]
    data: Optional[Any]
    created_at: datetime

    class Config:
        orm_mode = True
