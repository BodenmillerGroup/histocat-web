from typing import Dict, Optional

from pydantic import BaseModel


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
    class Config:
        orm_mode = True
