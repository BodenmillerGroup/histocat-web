import logging

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY

from histocat.db.base import Base

logger = logging.getLogger(__name__)


class ShareModel(Base):
    """Share."""

    __tablename__ = "share"
    __table_args__ = (sa.UniqueConstraint("user_id", "experiment_id", name="uq_share_user_id_experiment_id"),)

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    experiment_id = sa.Column(sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), index=True)
    user_id = sa.Column(sa.Integer(), sa.ForeignKey("user.id", ondelete="CASCADE"), index=True)
    permissions = sa.Column(ARRAY(sa.String(64)))
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    user = sa.orm.relationship("UserModel")
    experiment = sa.orm.relationship("ExperimentModel")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, user_id={self.user_id}, experiment_id={self.experiment_id})>"
