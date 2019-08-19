import logging
from datetime import datetime

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship

from app.db.base import Base

logger = logging.getLogger(__name__)


class Channel(Base):
    """
    Channel
    """

    __tablename__ = "channel"

    id: int = sa.Column('id', sa.Integer(), primary_key=True, index=True)
    acquisition_id: int = sa.Column('acquisition_id', sa.Integer(), sa.ForeignKey("acquisition.id", ondelete="CASCADE"),
                                    index=True)
    metaname: str = sa.Column('metaname', sa.String(4096))
    original_id: int = sa.Column('original_id', sa.Integer(), index=True)
    metal: str = sa.Column('metal', sa.String(), index=True)
    label: str = sa.Column('label', sa.String(), index=True)
    mass: int = sa.Column('mass', sa.Integer())
    max_intensity: float = sa.Column('max_intensity', sa.Float())
    min_intensity: float = sa.Column('min_intensity', sa.Float())
    meta: dict = sa.Column('meta', JSONB())
    location: str = sa.Column('location', sa.String(4096))
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    acquisition = relationship("Acquisition", back_populates="channels")

    def __repr__(self):
        return f"<Channel(id={self.id}, metal={self.metal}, label={self.label})>"
