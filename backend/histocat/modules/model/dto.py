from datetime import datetime
from typing import Optional, Dict

from pydantic import BaseModel


class ModelCreateDto(BaseModel):
    name: str
    description: Optional[str]


class ModelUpdateDto(BaseModel):
    name: str
    description: Optional[str]


class ModelDto(BaseModel):
    id: int
    group_id: int
    name: Optional[str]
    description: Optional[str]
    location: Optional[str]
    meta: Optional[Dict[str, Optional[str]]]
    created_at: datetime

    class Config:
        orm_mode = True
