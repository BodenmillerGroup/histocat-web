import logging
from datetime import datetime

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship

from app.db.base import Base

logger = logging.getLogger(__name__)


class ROIPoint(Base):
    """
    ROIPoint
    """

    __tablename__ = "roi_point"

    id: int = sa.Column(sa.Integer(), primary_key=True, index=True)
    roi_id: int = sa.Column(sa.Integer(), sa.ForeignKey("roi.id", ondelete="CASCADE"), index=True)
    metaname: str = sa.Column('metaname', sa.String(4096))
    original_id: int = sa.Column('original_id', sa.Integer(), index=True)
    order_number: int = sa.Column('order_number', sa.Integer())
    slide_x_pos_um: float = sa.Column('slide_x_pos_um', sa.Float())
    slide_y_pos_um: float = sa.Column('slide_y_pos_um', sa.Float())
    panorama_pixel_x_pos: int = sa.Column('panorama_pixel_x_pos', sa.Integer())
    panorama_pixel_y_pos: int = sa.Column('panorama_pixel_y_pos', sa.Integer())
    meta: dict = sa.Column('meta', JSONB())
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    roi = relationship("ROI", back_populates="roi_points")

    def __repr__(self):
        return f"<ROIPoint(id={self.id}, roi_id={self.roi_id}, order_number={self.order_number})>"
