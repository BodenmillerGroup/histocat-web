from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel


# Properties to receive via API on creation
class ShareCreateModel(BaseModel):
    experiment_id: int
    user_ids: List[int]
    permissions: Optional[List[str]]


# Shared properties
class ShareModel(BaseModel):
    id: int
    experiment_id: int
    user_id: int
    permissions: Optional[List[str]]
    created_at: datetime

    class Config:
        orm_mode = True
