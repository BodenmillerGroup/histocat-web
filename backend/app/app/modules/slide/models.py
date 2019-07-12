from datetime import datetime
from typing import List, Optional, Dict

from pydantic import BaseModel

from app.modules.panorama.models import PanoramaDatasetModel


# Properties to receive via API on creation
class SlideCreateModel(BaseModel):
    experiment_id: int
    metaname: str
    original_id: int
    uid: Optional[str]
    original_metadata: str
    meta: Dict[str, Optional[str]]


# Shared properties
class SlideModel(BaseModel):
    id: int
    experiment_id: int
    metaname: str
    original_id: int
    uid: Optional[str]
    original_metadata: str
    meta: Dict[str, Optional[str]]
    location: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


# Full slide dataset
class SlideDatasetModel(SlideModel):
    panoramas: List[PanoramaDatasetModel]

    class Config:
        orm_mode = True
