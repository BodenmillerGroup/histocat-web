from datetime import datetime
from typing import Optional

from pydantic import BaseModel, EmailStr


class UserCreateDto(BaseModel):
    email: EmailStr
    name: Optional[str]
    password: str
    is_active: Optional[bool] = True
    is_admin: Optional[bool] = False


class UserUpdateDto(BaseModel):
    email: EmailStr
    password: Optional[str]
    is_active: Optional[bool]
    is_admin: Optional[bool]
    name: Optional[str]


class UserDto(BaseModel):
    id: int
    email: EmailStr
    is_active: bool
    is_admin: bool
    name: Optional[str]
    created_at: datetime
    updated_at: datetime

    class Config:
        orm_mode = True
