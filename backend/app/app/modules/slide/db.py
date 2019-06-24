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

#: Format string for slide locations
SLIDE_LOCATION_FORMAT = "slide_{id}"


@remove_location_upon_delete
class Slide(Base):
    """
    Experiment slide
    """

    __tablename__ = "slide"
    __table_args__ = (UniqueConstraint("name"),)

    id: int = Column(Integer, primary_key=True, index=True)
    #: ID of parent experiment
    experiment_id: int = Column(
        Integer,
        ForeignKey("experiment.id", ondelete="CASCADE"),
        index=True,
    )
    #: location of the slide, e.g. absolute path to a directory on disk
    location: str = Column("location", String(4096))
    #: name given by user
    name: str = Column(String, index=True)
    #: original slide filename
    filename: str = Column(String(4096))
    #: slide width in μm
    width_um: int = Column(Integer)
    #: slide height in μm
    height_um: int = Column(Integer)
    #: description provided by user
    description: str = Column(Text)
    meta: dict = Column(JSONB)
    created_at: datetime = Column(DateTime, default=now(), nullable=False)

    experiment = relationship("Experiment", back_populates="slides")
    acquisitions = relationship("Acquisition", back_populates="slide", cascade="all, delete, delete-orphan")

    def __init__(
        self,
        experiment_id: int,
        name: str,
        filename: str,
        width_um: int,
        height_um: int,
        description: str = "",
        meta: dict = None,
    ):
        self.experiment_id = experiment_id
        self.name = name
        self.filename = filename
        self.width_um = width_um
        self.height_um = height_um
        self.description = description
        self.meta = meta

    @autocreate_directory_property
    def acquisitions_location(self) -> str:
        """
        Location where acquisitions are stored
        """
        return os.path.join(self.location, "acquisitions")

    def json(self):
        return {
            "id": self.id,
            "experiment_id": self.experiment_id,
            "name": self.name,
            "description": self.description,
            "location": self.location,
            "meta": self.meta,
            "created_at": self.created_at,
            "acquisitions": self.acquisitions,
        }

    def __repr__(self):
        return f"<Slide(id={self.id}, name={self.name})>"
