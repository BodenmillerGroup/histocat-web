from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel

from app.modules.slide.models import SlideDatasetModel


# Properties to receive via API on creation
class ExperimentCreateModel(BaseModel):
    name: str
    description: Optional[str] = None
    meta: Optional[dict] = None


# Properties to receive via API on update
class ExperimentUpdateModel(BaseModel):
    name: str
    description: Optional[str] = None
    meta: Optional[dict] = None


# Shared properties
class ExperimentBaseModel(BaseModel):
    id: int
    name: Optional[str]
    description: Optional[str]
    root_directory: str
    meta: Optional[dict]
    created_at: datetime


# Additional properties to return via API
class ExperimentModel(ExperimentBaseModel):
    location: Optional[str]


# Full experiment dataset
class ExperimentDatasetModel(ExperimentModel):
    slides: Optional[List[SlideDatasetModel]]
