import logging

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from histocat.core.utils import remove_location_upon_delete
from histocat.db.base import Base

logger = logging.getLogger(__name__)


@remove_location_upon_delete
class ResultModel(Base):
    """Result."""

    __tablename__ = "result"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    dataset_id = sa.Column(sa.Integer(), sa.ForeignKey("dataset.id", ondelete="CASCADE"), index=True, nullable=False)
    status = sa.Column(sa.String(64), nullable=False, index=True)
    name = sa.Column(sa.String())
    description = sa.Column(sa.String())
    pipeline = sa.Column(JSONB())
    input = sa.Column(JSONB())
    output = sa.Column(JSONB())
    location = sa.Column(sa.String(4096))
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    dataset = sa.orm.relationship("DatasetModel", back_populates="results")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
