from typing import Optional

from datetime import datetime
from pydantic import BaseModel


# Shared properties
class AcquisitionBaseModel(BaseModel):
    id: int = None
    slide_id: int = None
    name: Optional[str] = None
    description: Optional[str] = None
    location: Optional[str] = None
    meta: Optional[dict] = None


# Properties to receive via API on creation
class AcquisitionInCreateModel(BaseModel):
    name: str
    slide_id: int
    description: Optional[str] = None
    meta: Optional[dict] = None


# Properties to receive via API on update
class AcquisitionInUpdateModel(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    meta: Optional[dict] = None


# Additional properties to return via API
class AcquisitionModel(AcquisitionBaseModel):
    pass


# Additional properties stored in DB
class AcquisitionInDBModel(AcquisitionBaseModel):
    created_at: datetime
