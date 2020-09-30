from datetime import datetime
from typing import Any, Optional

from pydantic import BaseModel


class PipelineCreateDto(BaseModel):
    project_id: int
    name: Optional[str]
    description: Optional[str]
    steps: Optional[Any]


class PipelineUpdateDto(BaseModel):
    name: Optional[str]
    description: Optional[str]
    steps: Optional[Any]


class PipelineDto(BaseModel):
    id: int
    project_id: int
    name: Optional[str]
    description: Optional[str]
    steps: Optional[Any]
    created_at: datetime

    class Config:
        orm_mode = True
