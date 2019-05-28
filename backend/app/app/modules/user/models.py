from datetime import datetime
from typing import Optional

from pydantic import BaseModel, EmailStr, SecretStr


# Shared properties
class UserBaseModel(BaseModel):
    id: int = None
    email: Optional[EmailStr] = None
    is_active: Optional[bool] = True
    is_superuser: Optional[bool] = False
    full_name: Optional[str] = None


# Properties to receive via API on creation
class UserCreateModel(UserBaseModel):
    email: EmailStr
    password: SecretStr


# Properties to receive via API on update
class UserUpdateModel(UserBaseModel):
    password: Optional[SecretStr] = None


# Additional properties to return via API
class UserModel(UserBaseModel):
    created_at: datetime


# Additional properties stored in DB
class UserInDBModel(UserBaseModel):
    hashed_password: SecretStr
    created_at: datetime
