from datetime import datetime
from typing import Optional

from pydantic import BaseModel


class DatasetCreateModel(BaseModel):
    experiment_id: int
    user_id: int
    status: str
    name: Optional[str]
    description: Optional[str]
    artifacts: Optional[dict]
    meta: Optional[dict]
    errors: Optional[dict]


class DatasetUpdateModel(BaseModel):
    status: str
    name: Optional[str]
    description: Optional[str]
    artifacts: Optional[dict]
    meta: Optional[dict]
    errors: Optional[dict]


class DatasetModel(BaseModel):
    id: int
    experiment_id: int
    user_id: int
    uid: str
    status: str
    name: Optional[str]
    description: Optional[str]
    artifacts: Optional[dict]
    meta: Optional[dict]
    errors: Optional[dict]
    location: Optional[str]
    created_at: datetime
    updated_at: datetime

    class Config:
        orm_mode = True
