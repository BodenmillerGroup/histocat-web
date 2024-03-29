import os

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY

from histocat.core.base import Base
from histocat.core.utils import (
    autocreate_directory_property,
    remove_location_upon_delete,
)


@remove_location_upon_delete
class GroupModel(Base):
    """A *group* is the main organizational unit of `histoCAT`."""

    __tablename__ = "group"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    name = sa.Column(sa.String(), index=True, unique=True, nullable=False)
    description = sa.Column(sa.String())
    url = sa.Column(sa.String())
    is_open = sa.Column(sa.Boolean(), default=False, nullable=False)
    tags = sa.Column(ARRAY(sa.String(64)))
    location = sa.Column(sa.String(4096))
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    members = sa.orm.relationship("MemberModel", back_populates="group", cascade="all, delete, delete-orphan")
    projects = sa.orm.relationship("ProjectModel", back_populates="group", cascade="all, delete, delete-orphan")

    @autocreate_directory_property
    def projects_location(self) -> str:
        """Location where projects data are stored."""
        return os.path.join(self.location, "projects")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
