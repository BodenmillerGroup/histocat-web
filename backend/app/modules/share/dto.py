from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel


class ShareCreateDto(BaseModel):
    experiment_id: int
    user_ids: List[int]
    permissions: Optional[List[str]]


class ShareDto(BaseModel):
    id: int
    experiment_id: int
    user_id: int
    permissions: Optional[List[str]]
    created_at: datetime

    class Config:
        orm_mode = True
