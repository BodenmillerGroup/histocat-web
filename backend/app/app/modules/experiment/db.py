import logging
import os
from datetime import datetime
from typing import List

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY, JSONB
from sqlalchemy.orm import relationship

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

    id: int = sa.Column(sa.Integer(), primary_key=True, index=True)
    user_id: int = sa.Column(sa.Integer(), sa.ForeignKey("user.id", ondelete="CASCADE"), index=True)
    name: str = sa.Column(sa.String(), index=True)
    description: str = sa.Column(sa.Text())
    meta: dict = sa.Column(JSONB())
    tags: List[str] = sa.Column(ARRAY(sa.String(64)))
    location: str = sa.Column("location", sa.String(4096))
    created_at: datetime = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    user = relationship("User", back_populates="experiments")
    slides = relationship("Slide", back_populates="experiment", cascade="all, delete, delete-orphan")
    datasets = relationship("Dataset", back_populates="experiment", cascade="all, delete, delete-orphan")

    @autocreate_directory_property
    def slides_location(self) -> str:
        """
        Location where slides data are stored
        """
        return os.path.join(self.location, "slides")

    @autocreate_directory_property
    def datasets_location(self) -> str:
        """
        Location where datasets are stored
        """
        return os.path.join(self.location, "datasets")

    def __repr__(self):
        return f"<Experiment(id={self.id}, name={self.name})>"
