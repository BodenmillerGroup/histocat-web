from typing import Optional

from datetime import datetime
from pydantic import BaseModel


# Shared properties
class ChannelBaseModel(BaseModel):
    id: int = None
    acquisition_id: int = None
    name: Optional[str] = None
    metal: Optional[str] = None
    location: Optional[str] = None
    meta: Optional[dict] = None


# Properties to receive via API on creation
class ChannelInCreateModel(BaseModel):
    name: str
    metal: str
    acquisition_id: int
    meta: Optional[dict] = None


# Properties to receive via API on update
class ChannelInUpdateModel(BaseModel):
    name: Optional[str] = None
    meta: Optional[dict] = None


# Additional properties to return via API
class ChannelModel(ChannelBaseModel):
    pass


# Additional properties stored in DB
class ChannelInDBModel(ChannelBaseModel):
    created_at: datetime
