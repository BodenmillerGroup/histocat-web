import logging
import os
from datetime import datetime

from sqlalchemy import Column, ForeignKey, Integer, String, Text, UniqueConstraint, DateTime
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship
from sqlalchemy.sql.functions import now

from app.core.utils import autocreate_directory_property, remove_location_upon_delete
from app.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for acquisition locations
ACQUISITION_LOCATION_FORMAT = "acquisition_{id}"


@remove_location_upon_delete
class Acquisition(Base):
    """
    An *acquisition* contains all files belonging to one microscope image acquisition process.
    """

    __tablename__ = "acquisition"
    __table_args__ = (UniqueConstraint("name", "slide_id"),)

    id: int = Column(Integer, primary_key=True, index=True)
    #: ID of the parent slide
    slide_id: int = Column(
        Integer,
        ForeignKey("slide.id", onupdate="CASCADE", ondelete="CASCADE"),
        index=True,
    )
    #: location of the acquisition, e.g. absolute path to a directory on disk
    location: str = Column("location", String(4096))
    #: name
    name: str = Column(String, index=True)
    #: acquisition width in pixels
    width: int = Column(Integer)
    #: acquisition height in pixels
    height: int = Column(Integer)
    #: description
    description: str = Column(Text)
    meta: dict = Column(JSONB)
    created_at: datetime = Column(DateTime, default=now(), nullable=False)

    slide = relationship("Slide", back_populates="acquisitions")
    channels = relationship("Channel", back_populates="acquisition")

    def __init__(
        self,
        slide_id: int,
        name: str,
        width: int,
        height: int,
        description: str = "",
        meta: dict = None,
    ):
        # TODO: ensure that name is unique within slide
        self.slide_id = slide_id
        self.name = name
        self.width = width
        self.height = height
        self.description = description
        self.meta = meta

    @autocreate_directory_property
    def channels_location(self) -> str:
        """str: location where channels files are stored"""
        return os.path.join(self.location, "channels")

    def to_dict(self) -> dict:
        """Returns attributes as key-value pairs.

        Returns
        -------
        dict
        """
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "status": self.status,
        }

    def json(self) -> dict:
        return {
            "id": self.id,
            "slide_id": self.slide_id,
            "name": self.name,
            "description": self.description,
            "location": self.location,
            "meta": self.meta,
            "created_at": self.created_at,
            "channels": self.channels,
        }

    def __repr__(self):
        return f"<Acquisition(id={self.id}, name={self.name})>"
