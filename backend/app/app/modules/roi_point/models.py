from datetime import datetime
from typing import Dict, Optional

from pydantic import BaseModel


# Properties to receive via API on creation
class ROIPointCreateModel(BaseModel):
    roi_id: int
    origin_id: int
    meta: Dict[str, Optional[str]]


# Shared properties
class ROIPointModel(BaseModel):
    id: int
    roi_id: int
    origin_id: int
    meta: Dict[str, Optional[str]]
    created_at: datetime

    class Config:
        orm_mode = True
