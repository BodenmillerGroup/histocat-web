import logging
import os

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY

from histocat.core.utils import (
    autocreate_directory_property,
    remove_location_upon_delete,
)
from histocat.db.base import Base

logger = logging.getLogger(__name__)


@remove_location_upon_delete
class ProjectModel(Base):
    """An *project* is the secondary organizational unit of `histoCAT`."""

    __tablename__ = "project"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    group_id = sa.Column(sa.Integer(), sa.ForeignKey("group.id", ondelete="CASCADE"), index=True, nullable=False)
    member_id = sa.Column(sa.Integer(), sa.ForeignKey("member.id", ondelete="CASCADE"), index=True, nullable=False)
    name = sa.Column(sa.String(), index=True)
    description = sa.Column(sa.Text())
    tags = sa.Column(ARRAY(sa.String(64)))
    location = sa.Column(sa.String(4096))
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    group = sa.orm.relationship("GroupModel", back_populates="projects")
    member = sa.orm.relationship("MemberModel", back_populates="projects")
    slides = sa.orm.relationship("SlideModel", back_populates="project", cascade="all, delete, delete-orphan")
    datasets = sa.orm.relationship("DatasetModel", back_populates="project", cascade="all, delete, delete-orphan")
    presets = sa.orm.relationship("PresetModel", back_populates="project", cascade="all, delete, delete-orphan")

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
