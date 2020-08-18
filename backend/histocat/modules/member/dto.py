from datetime import datetime

from pydantic import BaseModel


class MemberCreateDto(BaseModel):
    group_id: int
    user_id: int
    role: int
    is_active: bool


class MemberUpdateDto(BaseModel):
    role: int
    is_active: bool


class MemberDto(BaseModel):
    id: int
    group_id: int
    user_id: int
    role: int
    is_active: bool
    created_at: datetime
    updated_at: datetime

    class Config:
        orm_mode = True
