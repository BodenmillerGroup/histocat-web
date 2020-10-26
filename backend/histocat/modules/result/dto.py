from datetime import datetime
from typing import Any, Optional

from pydantic import BaseModel


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
