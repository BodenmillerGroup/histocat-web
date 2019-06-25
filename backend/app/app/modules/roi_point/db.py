import logging
from datetime import datetime

import sqlalchemy as sa
from sqlalchemy.orm import relationship

from app.core.utils import remove_location_upon_delete
from app.db.base import Base

logger = logging.getLogger(__name__)


@remove_location_upon_delete
class ROIPoint(Base):
    """
    ROIPoint
    """

    __tablename__ = "roi_point"

    id: int = sa.Column(sa.Integer(), primary_key=True, index=True)
    roi_id: int = sa.Column(sa.Integer(), sa.ForeignKey("roi.id", ondelete="CASCADE"), index=True)
    order_number: int = sa.Column('order_number', sa.Integer())
    slide_x_pos_um: float = sa.Column('slide_x_pos_um', sa.Float())
    slide_y_pos_um: float = sa.Column('slide_y_pos_um', sa.Float())
    panorama_pixel_x_pos: int = sa.Column('panorama_pixel_x_pos', sa.Integer())
    panorama_pixel_y_pos: int = sa.Column('panorama_pixel_y_pos', sa.Integer())
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    roi = relationship("ROI", back_populates="roi_points")

    def __repr__(self):
        return f"<ROIPoint(id={self.id}, roi_id={self.roi_id}, order_number={self.order_number})>"
