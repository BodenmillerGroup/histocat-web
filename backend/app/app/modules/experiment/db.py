import logging
import os
from datetime import datetime

from sqlalchemy import Column, String, Text, Integer, ForeignKey, DateTime
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship
from sqlalchemy.sql.functions import now

from app.core.utils import autocreate_directory_property, remove_location_upon_delete
from app.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for experiment locations.
EXPERIMENT_LOCATION_FORMAT = "experiment_{id}"


@remove_location_upon_delete
class Experiment(Base):
    """
    An *experiment* is the main organizational unit of `HistoCAT`.
    """

    __tablename__ = "experiment"

    id: int = Column(Integer, primary_key=True, index=True)
    #: ID of the user
    user_id: int = Column(
        Integer,
        ForeignKey("user.id", onupdate="CASCADE", ondelete="CASCADE"),
        index=True,
    )
    #: location of the experiment, e.g. absolute path to a directory on disk
    location: str = Column("location", String(4096))
    #: experiment name
    name: str = Column(String, index=True)
    #: experiment description
    description: str = Column(Text)
    meta: dict = Column(JSONB)
    created_at: datetime = Column(DateTime, default=now(), nullable=False)

    user = relationship("User", back_populates="experiments")
    slides = relationship("Slide", back_populates="experiment")

    def __init__(self, name: str, user_id: int, description: str = "", meta: dict = None):
        self.name = name
        self.user_id = user_id
        self.description = description
        self.meta = meta

    @autocreate_directory_property
    def slides_location(self) -> str:
        """
        Location where slides data are stored
        """
        return os.path.join(self.location, "slides")

    def json(self) -> dict:
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "location": self.location,
            "meta": self.meta,
            "user_id": self.user_id,
            "created_at": self.created_at,
            "slides": self.slides,
        }

    def __repr__(self):
        return f"<Experiment(id={self.id}, name={self.name})>"
