from datetime import datetime
from typing import Optional, Sequence

from pydantic import BaseModel


class GateCreateDto(BaseModel):
    dataset_id: int
    name: Optional[str]
    description: Optional[str]
    indices: Sequence[int]
    cell_ids: Sequence[int]


class GateDto(BaseModel):
    id: int
    dataset_id: int
    name: Optional[str]
    description: Optional[str]
    indices: Sequence[int]
    cell_ids: Sequence[int]
    created_at: datetime

    class Config:
        orm_mode = True
