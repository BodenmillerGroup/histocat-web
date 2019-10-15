from datetime import datetime
from typing import List, Optional, Dict

from pydantic import BaseModel

from app.modules.panorama.models import PanoramaDatasetModel


class SlideCreateModel(BaseModel):
    experiment_id: int
    name: str
    origin_id: int
    xml_meta: str
    meta: Dict[str, Optional[str]]


class SlideModel(BaseModel):
    id: int
    experiment_id: int
    name: str
    origin_id: int
    meta: Dict[str, Optional[str]]
    created_at: datetime

    class Config:
        orm_mode = True


class SlideDatasetModel(SlideModel):
    """
    Full slide dataset
    """

    panoramas: List[PanoramaDatasetModel]

    class Config:
        orm_mode = True
