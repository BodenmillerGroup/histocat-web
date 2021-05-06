from datetime import datetime
from typing import Any, Optional

from pydantic import BaseModel


class GateCreateDto(BaseModel):
    dataset_id: int
    name: Optional[str]
    description: Optional[str]
    cell_classes: Optional[Any]
    annotations: Optional[Any]


class GateUpdateDto(BaseModel):
    name: Optional[str]
    description: Optional[str]


class GateDto(BaseModel):
    id: int
    dataset_id: int
    name: Optional[str]
    description: Optional[str]
    cell_classes: Optional[Any]
    annotations: Optional[Any]
    created_at: datetime

    class Config:
        orm_mode = True
