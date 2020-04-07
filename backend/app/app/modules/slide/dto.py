from datetime import datetime
from typing import Dict, Optional, Sequence

from pydantic import BaseModel

from app.modules.acquisition.dto import AcquisitionDatasetDto
from app.modules.panorama.dto import PanoramaDatasetDto


class SlideCreateDto(BaseModel):
    experiment_id: int
    name: str
    origin_id: int
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
    panoramas: Sequence[PanoramaDatasetDto]
    acquisitions: Sequence[AcquisitionDatasetDto]

    class Config:
        orm_mode = True
