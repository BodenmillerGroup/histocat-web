from datetime import datetime
from typing import Optional, Any

from pydantic import BaseModel


class ResultCreateDto(BaseModel):
    dataset_id: int
    status: str
    name: Optional[str]
    description: Optional[str]
    pipeline: Any
    input: Any


class ResultUpdateDto(BaseModel):
    dataset_id: int
    status: Optional[str]
    name: Optional[str]
    description: Optional[str]


class ResultDto(BaseModel):
    id: int
    dataset_id: int
    status: str
    name: Optional[str]
    description: Optional[str]
    pipeline: Any
    input: Any
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True
