from datetime import datetime
from typing import Optional

from pydantic import BaseModel


# Properties to receive via API on creation
class ROICreateModel(BaseModel):
    panorama_id: int
    roi_type: Optional[str]


# Shared properties
class ROIModel(BaseModel):
    id: int
    panorama_id: int
    roi_type: Optional[str]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True
