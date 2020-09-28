from datetime import datetime
from typing import Optional

from pydantic import BaseModel


class ResultCreateDto(BaseModel):
    dataset_id: int
    parent_id: Optional[int]
    type: str
    status: str
    name: Optional[str]
    description: Optional[str]
    params: dict


class ResultUpdateDto(BaseModel):
    dataset_id: int
    status: Optional[str]
    name: Optional[str]
    description: Optional[str]


class ResultDto(BaseModel):
    id: int
    dataset_id: int
    parent_id: Optional[int]
    type: str
    status: str
    name: Optional[str]
    description: Optional[str]
    params: dict
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True
