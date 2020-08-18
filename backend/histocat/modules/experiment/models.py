import logging
import os

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY, JSONB

from histocat.core.utils import (
    autocreate_directory_property,
    remove_location_upon_delete,
)
from histocat.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for experiment locations.
EXPERIMENT_LOCATION_FORMAT = "experiment_{id}"


@remove_location_upon_delete
class ExperimentModel(Base):
    """An *experiment* is the main organizational unit of `HistoCAT`."""

    __tablename__ = "experiment"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    group_id = sa.Column(sa.Integer(), sa.ForeignKey("group.id", ondelete="CASCADE"), index=True)
    name = sa.Column(sa.String(), index=True)
    description = sa.Column(sa.Text())
    meta = sa.Column(JSONB())
    tags = sa.Column(ARRAY(sa.String(64)))
    is_public = sa.Column(sa.Boolean(), nullable=False, default=False)
    location = sa.Column(sa.String(4096))
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    group = sa.orm.relationship("GroupModel", back_populates="experiments")
    slides = sa.orm.relationship("SlideModel", back_populates="experiment", cascade="all, delete, delete-orphan")
    datasets = sa.orm.relationship("DatasetModel", back_populates="experiment", cascade="all, delete, delete-orphan")
    presets = sa.orm.relationship("PresetModel", back_populates="experiment", cascade="all, delete, delete-orphan")

    @autocreate_directory_property
    def slides_location(self) -> str:
        """Location where slides data are stored."""
        return os.path.join(self.location, "slides")

    @autocreate_directory_property
    def datasets_location(self) -> str:
        """Location where datasets are stored."""
        return os.path.join(self.location, "datasets")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
