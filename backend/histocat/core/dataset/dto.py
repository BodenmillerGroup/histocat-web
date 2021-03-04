from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel


class DatasetCreateDto(BaseModel):
    project_id: int
    status: str
    name: Optional[str]
    description: Optional[str]
    origin: str
    acquisition_ids: Optional[List[int]]
    channels: Optional[List[str]]
    meta: Optional[dict]


class DatasetUpdateDto(BaseModel):
    status: Optional[str]
    name: Optional[str]
    description: Optional[str]
    acquisition_ids: Optional[List[int]]
    channels: Optional[List[str]]
    meta: Optional[dict]


class DatasetDto(BaseModel):
    id: int
    project_id: int
    uid: str
    status: str
    name: Optional[str]
    description: Optional[str]
    origin: str
    acquisition_ids: Optional[List[int]]
    channels: Optional[List[str]]
    meta: Optional[dict]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True
