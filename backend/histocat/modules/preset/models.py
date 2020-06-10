import logging

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from histocat.db.base import Base

logger = logging.getLogger(__name__)


class PresetModel(Base):
    """User custom experiment preset."""

    __tablename__ = "preset"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    experiment_id = sa.Column(sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), index=True)
    name = sa.Column(sa.String())
    description = sa.Column(sa.String())
    data = sa.Column(JSONB())
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    experiment = sa.orm.relationship("ExperimentModel", back_populates="presets")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
