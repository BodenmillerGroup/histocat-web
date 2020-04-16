from datetime import datetime
from typing import Optional

from pydantic import BaseModel


class DatasetCreateDto(BaseModel):
    experiment_id: int
    user_id: int
    status: str
    name: Optional[str]
    description: Optional[str]
    input: Optional[dict]
    meta: Optional[dict]


class DatasetUpdateDto(BaseModel):
    status: str
    name: Optional[str]
    description: Optional[str]
    input: Optional[dict]
    output: Optional[dict]
    meta: Optional[dict]


class DatasetDto(BaseModel):
    id: int
    experiment_id: int
    user_id: int
    uid: str
    status: str
    name: Optional[str]
    description: Optional[str]
    input: Optional[dict]
    output: Optional[dict]
    meta: Optional[dict]
    location: Optional[str]
    created_at: datetime
    updated_at: datetime

    class Config:
        orm_mode = True