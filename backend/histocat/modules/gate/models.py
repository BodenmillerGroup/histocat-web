import logging

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY

from histocat.db.base import Base

logger = logging.getLogger(__name__)


class GateModel(Base):
    """Cell gate."""

    __tablename__ = "gate"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    dataset_id = sa.Column(sa.Integer(), sa.ForeignKey("dataset.id", ondelete="CASCADE"), index=True)
    name = sa.Column(sa.String())
    description = sa.Column(sa.String())
    acquisition_ids = sa.Column("acquisition_ids", ARRAY(sa.Integer()))
    indices = sa.Column("indices", ARRAY(sa.Integer()))
    cell_ids = sa.Column("cell_ids", ARRAY(sa.Integer()))
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    dataset = sa.orm.relationship("DatasetModel", back_populates="gates")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
