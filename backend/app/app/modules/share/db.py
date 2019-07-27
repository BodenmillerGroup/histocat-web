import logging
from datetime import datetime
from typing import List

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import relationship

from app.db.base import Base

logger = logging.getLogger(__name__)


class Share(Base):
    """
    Share
    """

    __tablename__ = "share"
    __table_args__ = (
        sa.UniqueConstraint('user_id', 'experiment_id', name='uq_share_user_id_experiment_id'),
    )

    id: int = sa.Column(sa.Integer(), primary_key=True, index=True)
    experiment_id: int = sa.Column(sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), index=True)
    user_id: int = sa.Column(sa.Integer(), sa.ForeignKey("user.id", ondelete="CASCADE"), index=True)
    permissions: List[str] = sa.Column(ARRAY(sa.String(64)))
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    user = relationship("User")
    experiment = relationship("Experiment")

    def __repr__(self):
        return f"<Share(id={self.id}, user_id={self.user_id}, experiment_id={self.experiment_id})>"
