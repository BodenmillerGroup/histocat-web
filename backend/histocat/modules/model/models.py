import logging

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from histocat.core.utils import remove_location_upon_delete
from histocat.db.base import Base

logger = logging.getLogger(__name__)


@remove_location_upon_delete
class ModelModel(Base):
    """Segmentation Keras model."""

    __tablename__ = "model"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    group_id = sa.Column(sa.Integer(), sa.ForeignKey("group.id", ondelete="CASCADE"), index=True, nullable=False)
    name = sa.Column(sa.String(), index=True)
    description = sa.Column(sa.Text())
    location = sa.Column(sa.String(4096))
    meta = sa.Column(JSONB())
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    group = sa.orm.relationship("GroupModel", back_populates="models")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
