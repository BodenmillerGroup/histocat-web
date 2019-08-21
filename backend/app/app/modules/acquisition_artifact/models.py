from datetime import datetime
from typing import Optional, Dict

from pydantic import BaseModel


# Properties to receive via API on creation
class AcquisitionArtifactCreateModel(BaseModel):
    acquisition_id: int
    type: str
    description: Optional[str]
    location: Optional[str]
    meta: Dict[str, Optional[str]]


# Shared properties
class AcquisitionArtifactModel(BaseModel):
    id: int
    acquisition_id: int
    type: str
    description: Optional[str]
    location: Optional[str]
    meta: Dict[str, Optional[str]]
    created_at: datetime

    class Config:
        orm_mode = True
