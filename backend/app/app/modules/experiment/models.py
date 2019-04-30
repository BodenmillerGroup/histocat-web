from datetime import datetime
from typing import Optional, List

from pydantic import BaseModel

from app.modules.slide.models import SlideDatasetModel


# Shared properties
class ExperimentBaseModel(BaseModel):
    id: int = None
    name: Optional[str] = None
    description: Optional[str] = None
    root_directory: Optional[str] = None
    meta: Optional[dict] = None
    created_at: datetime


# Properties to receive via API on creation
class ExperimentCreateModel(BaseModel):
    name: str
    description: Optional[str] = None
    meta: Optional[dict] = None


# Properties to receive via API on update
class ExperimentUpdateModel(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    meta: Optional[dict] = None


# Additional properties to return via API
class ExperimentModel(ExperimentBaseModel):
    location: Optional[str] = None


# Additional properties stored in DB
class ExperimentInDBModel(ExperimentBaseModel):
    pass


# Full experiment dataset
class ExperimentDatasetModel(ExperimentModel):
    slides: Optional[List[SlideDatasetModel]]
