from typing import Optional

from datetime import datetime
from pydantic import BaseModel


# Shared properties
class ExperimentBaseModel(BaseModel):
    id: int = None
    name: Optional[str] = None
    description: Optional[str] = None
    root_directory: Optional[str] = None
    location: Optional[str] = None


# Properties to receive via API on creation
class ExperimentInCreateModel(ExperimentBaseModel):
    name: str
    description: Optional[str] = None
    root_directory: str


# Properties to receive via API on update
class ExperimentInUpdateModel(ExperimentBaseModel):
    name: Optional[str] = None
    description: Optional[str] = None


# Additional properties to return via API
class ExperimentModel(ExperimentBaseModel):
    pass


# Additional properties stored in DB
class ExperimentInDBModel(ExperimentBaseModel):
    meta: dict
    created_at: datetime
