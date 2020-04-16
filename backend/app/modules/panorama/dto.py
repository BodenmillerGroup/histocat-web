from typing import Dict, Optional

from pydantic import BaseModel


class PanoramaCreateDto(BaseModel):
    slide_id: int
    origin_id: int
    image_type: Optional[str]
    description: Optional[str]
    start_position_x: Optional[float]
    start_position_y: Optional[float]
    width: Optional[float]
    height: Optional[float]
    rotation_angle: Optional[float]
    location: Optional[str]
    meta: Dict[str, Optional[str]]


class PanoramaDto(BaseModel):
    id: int
    slide_id: int
    origin_id: int
    image_type: Optional[str]
    description: Optional[str]
    start_position_x: Optional[float]
    start_position_y: Optional[float]
    width: Optional[float]
    height: Optional[float]
    rotation_angle: Optional[float]
    meta: Dict[str, Optional[str]]

    class Config:
        orm_mode = True


class PanoramaDatasetDto(PanoramaDto):
    """Full panorama dataset."""

    class Config:
        orm_mode = True
