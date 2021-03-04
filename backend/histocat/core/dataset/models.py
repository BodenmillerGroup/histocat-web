import os

import sqlalchemy as sa
from sqlalchemy import text
from sqlalchemy.dialects.postgresql import ARRAY, JSONB, UUID

from histocat.core.base import Base
from histocat.core.constants import ANNDATA_FILE_EXTENSION
from histocat.core.utils import (
    autocreate_directory_property,
    remove_location_upon_delete,
)

CELL_FILENAME = f"cell{ANNDATA_FILE_EXTENSION}"


@remove_location_upon_delete
class DatasetModel(Base):
    """Dataset."""

    __tablename__ = "dataset"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    project_id = sa.Column(sa.Integer(), sa.ForeignKey("project.id", ondelete="CASCADE"), index=True, nullable=False)
    uid = sa.Column(UUID(), server_default=text("uuid_generate_v4()"), index=True, nullable=False)
    status = sa.Column(sa.String(64), index=True, nullable=False)
    name = sa.Column(sa.String())
    description = sa.Column(sa.String())
    origin = sa.Column(sa.String())
    acquisition_ids = sa.Column("acquisition_ids", ARRAY(sa.Integer()))
    channels = sa.Column("channels", ARRAY(sa.String()))
    meta = sa.Column(JSONB())
    location = sa.Column(sa.String(4096))
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    project = sa.orm.relationship("ProjectModel", back_populates="datasets")
    gates = sa.orm.relationship("GateModel", back_populates="dataset", cascade="all, delete, delete-orphan")
    results = sa.orm.relationship("ResultModel", back_populates="dataset", cascade="all, delete, delete-orphan")

    @autocreate_directory_property
    def results_location(self) -> str:
        """Location where results are stored."""
        return os.path.join(self.location, "results")

    def cell_file_location(self) -> str:
        """Location where single-cell data are stored."""
        return os.path.join(self.location, CELL_FILENAME)

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
