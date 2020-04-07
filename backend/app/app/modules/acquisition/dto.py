from typing import Dict, Optional, Sequence

from pydantic import BaseModel

from app.modules.channel.dto import ChannelDto


class AcquisitionCreateDto(BaseModel):
    slide_id: int
    origin_id: int
    location: str
    meta: Dict[str, Optional[str]]


class AcquisitionDto(BaseModel):
    id: int
    slide_id: int
    origin_id: int
    meta: Dict[str, Optional[str]]

    class Config:
        orm_mode = True


class AcquisitionDatasetDto(AcquisitionDto):
    channels: Optional[Sequence[ChannelDto]]

    class Config:
        orm_mode = True
