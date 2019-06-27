from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel

from app.modules.slide.models import SlideDatasetModel


# Properties to receive via API on creation
class ExperimentCreateModel(BaseModel):
    name: str
    description: Optional[str]
    meta: Optional[dict]
    tags: Optional[List[str]]


# Properties to receive via API on update
class ExperimentUpdateModel(BaseModel):
    name: str
    description: Optional[str]
    meta: Optional[dict]
    tags: Optional[List[str]]


# Shared properties
class ExperimentModel(BaseModel):
    id: int
    user_id: int
    name: Optional[str]
    description: Optional[str]
    meta: Optional[dict]
    tags: Optional[List[str]]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


# Full experiment dataset
class ExperimentDatasetModel(ExperimentModel):
    slides: List[SlideDatasetModel]

    class Config:
        orm_mode = True
