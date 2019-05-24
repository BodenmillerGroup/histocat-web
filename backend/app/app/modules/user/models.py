from datetime import datetime
from typing import Optional

from pydantic import BaseModel


# Shared properties
class UserBaseModel(BaseModel):
    id: int = None
    email: Optional[str] = None
    is_active: Optional[bool] = True
    is_superuser: Optional[bool] = False
    full_name: Optional[str] = None


# Properties to receive via API on creation
class UserCreateModel(UserBaseModel):
    email: str
    password: str


# Properties to receive via API on update
class UserUpdateModel(UserBaseModel):
    password: Optional[str] = None


# Additional properties to return via API
class UserModel(UserBaseModel):
    pass


# Additional properties stored in DB
class UserInDBModel(UserBaseModel):
    hashed_password: str
    created_at: datetime
