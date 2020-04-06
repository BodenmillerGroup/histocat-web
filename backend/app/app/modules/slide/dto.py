from datetime import datetime
from typing import Dict, List, Optional

from pydantic import BaseModel

from app.modules.panorama.dto import PanoramaDatasetDto


class SlideCreateDto(BaseModel):
    experiment_id: int
    name: str
    origin_id: int
    xml_meta: str
    meta: Dict[str, Optional[str]]


class SlideDto(BaseModel):
    id: int
    experiment_id: int
    name: str
    origin_id: int
    meta: Dict[str, Optional[str]]
    created_at: datetime

    class Config:
        orm_mode = True


class SlideDatasetDto(SlideDto):
    """Full slide dataset."""
    panoramas: List[PanoramaDatasetDto]

    class Config:
        orm_mode = True
