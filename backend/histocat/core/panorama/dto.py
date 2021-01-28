from typing import Dict, Optional

from pydantic import BaseModel


class PanoramaCreateDto(BaseModel):
    slide_id: int
    origin_id: int
    image_type: Optional[str]
    description: Optional[str]
    x1: Optional[float]
    y1: Optional[float]
    x2: Optional[float]
    y2: Optional[float]
    x3: Optional[float]
    y3: Optional[float]
    x4: Optional[float]
    y4: Optional[float]
    rotation_angle: Optional[float]
    location: Optional[str]
    meta: Dict[str, Optional[str]]


class PanoramaDto(BaseModel):
    id: int
    slide_id: int
    origin_id: int
    image_type: Optional[str]
    description: Optional[str]
    x1: Optional[float]
    y1: Optional[float]
    x2: Optional[float]
    y2: Optional[float]
    x3: Optional[float]
    y3: Optional[float]
    x4: Optional[float]
    y4: Optional[float]
    rotation_angle: Optional[float]
    meta: Dict[str, Optional[str]]

    class Config:
        orm_mode = True


class PanoramaDatasetDto(PanoramaDto):
    """Full panorama dataset."""

    class Config:
        orm_mode = True
