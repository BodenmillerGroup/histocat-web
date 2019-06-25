import logging
from datetime import datetime

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship

from app.core.utils import remove_location_upon_delete
from app.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for channel locations
CHANNEL_LOCATION_FORMAT = "channel_{id}"


@remove_location_upon_delete
class Channel(Base):
    """
    Channel
    """

    __tablename__ = "channel"

    id: int = sa.Column('id', sa.Integer(), primary_key=True, index=True)
    acquisition_id: int = sa.Column('acquisition_id', sa.Integer(), sa.ForeignKey("acquisition.id", ondelete="CASCADE"),
                                    index=True)
    channel_name: str = sa.Column('channel_name', sa.String(), index=True)
    channel_label: str = sa.Column('channel_label', sa.String(), index=True)
    order_number: int = sa.Column('order_number', sa.Integer())
    mass: int = sa.Column('mass', sa.Integer())
    max_intensity: int = sa.Column('max_intensity', sa.Integer())
    min_intensity: int = sa.Column('min_intensity', sa.Integer())
    meta: dict = sa.Column('meta', JSONB())
    location: str = sa.Column('location', sa.String(4096))
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    acquisition = relationship("Acquisition", back_populates="channels")

    def __repr__(self):
        return f"<Channel(id={self.id}, channel_name={self.channel_name}, channel_label={self.channel_label})>"
