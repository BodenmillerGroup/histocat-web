from datetime import datetime
from typing import Dict

from pydantic import BaseModel


# Properties to receive via API on creation
class ROIPointCreateModel(BaseModel):
    roi_id: int
    metaname: str
    original_id: int
    meta: Dict[str, str]


# Shared properties
class ROIPointModel(BaseModel):
    id: int
    roi_id: int
    metaname: str
    original_id: int
    meta: Dict[str, str]
    created_at: datetime

    class Config:
        orm_mode = True
