from datetime import datetime
from typing import Optional

from pydantic import BaseModel


# Properties to receive via API on creation
class DatasetCreateModel(BaseModel):
    experiment_id: int
    name: str
    description: Optional[str]
    input: Optional[dict]
    meta: Optional[dict]


# Shared properties
class DatasetModel(BaseModel):
    id: int
    experiment_id: int
    user_id: int
    uid: str
    name: str
    description: Optional[str]
    status: str
    errors: Optional[dict]
    input: Optional[dict]
    meta: Optional[dict]
    location: Optional[str]
    created_at: datetime
    updated_at: datetime

    class Config:
        orm_mode = True
