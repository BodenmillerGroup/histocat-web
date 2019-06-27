from datetime import datetime
from typing import Optional

from pydantic import BaseModel, EmailStr, SecretStr


# Properties to receive via API on creation
class UserCreateModel(BaseModel):
    email: EmailStr
    password: str
    is_active: Optional[bool] = True
    is_superuser: Optional[bool] = False
    full_name: Optional[str]


# Properties to receive via API on update
class UserUpdateModel(BaseModel):
    email: EmailStr
    password: Optional[str]
    is_active: Optional[bool]
    is_superuser: Optional[bool]
    full_name: Optional[str]


# Shared properties
class UserModel(BaseModel):
    id: int
    email: EmailStr
    is_active: bool
    is_superuser: bool
    full_name: Optional[str]
    created_at: datetime

    class Config:
        orm_mode = True


# Additional properties stored in DB
class UserDBModel(UserModel):
    hashed_password: SecretStr
