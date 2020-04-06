from typing import Dict, List, Optional

from pydantic import BaseModel

from app.modules.roi.models import ROIDatasetModel


class PanoramaCreateDto(BaseModel):
    slide_id: int
    origin_id: int
    meta: Dict[str, Optional[str]]


class PanoramaDto(BaseModel):
    id: int
    slide_id: int
    origin_id: int
    meta: Dict[str, Optional[str]]

    class Config:
        orm_mode = True


class PanoramaDatasetDto(PanoramaDto):
    """Full panorama dataset."""
    rois: List[ROIDatasetModel]

    class Config:
        orm_mode = True
