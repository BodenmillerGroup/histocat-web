from typing import Optional

from datetime import datetime
from pydantic import BaseModel


# Shared properties
class AcquisitionBaseModel(BaseModel):
    id: int = None
    slide_id: int = None
    name: Optional[str] = None
    width: Optional[int] = None
    height: Optional[int] = None
    description: Optional[str] = None
    meta: Optional[dict] = None


# Properties to receive via API on creation
class AcquisitionCreateModel(BaseModel):
    slide_id: int
    name: str
    width: int
    height: int
    description: Optional[str] = None
    meta: Optional[dict] = None


# Properties to receive via API on update
class AcquisitionUpdateModel(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    meta: Optional[dict] = None


# Additional properties to return via API
class AcquisitionModel(AcquisitionBaseModel):
    location: Optional[str] = None


# Additional properties stored in DB
class AcquisitionInDBModel(AcquisitionBaseModel):
    created_at: datetime
