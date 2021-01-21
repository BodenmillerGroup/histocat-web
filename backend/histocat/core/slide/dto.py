from datetime import datetime
from typing import Any, Dict, Optional, Sequence

from pydantic import BaseModel

from histocat.core.acquisition.dto import AcquisitionDatasetDto
from histocat.core.panorama.dto import PanoramaDatasetDto


class SlideCreateDto(BaseModel):
    project_id: int
    origin_id: int
    name: str
    width_um: Optional[int]
    height_um: Optional[int]
    has_slide_image: Optional[bool]
    meta: Dict[str, Optional[str]]
    session_meta: Optional[Dict[str, Any]]


class SlideDto(BaseModel):
    id: int
    project_id: int
    origin_id: int
    name: str
    width_um: Optional[int]
    height_um: Optional[int]
    has_slide_image: bool
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
