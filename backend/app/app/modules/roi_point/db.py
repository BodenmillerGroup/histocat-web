import logging
from datetime import datetime
from typing import Optional

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
    meta: dict = sa.Column('meta', JSONB())
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    roi = relationship("ROI", back_populates="roi_points")

    @property
    def OrderNumber(self) -> Optional[str]:
        return self.meta.get('OrderNumber')

    @property
    def SlideXPosUm(self) -> Optional[str]:
        return self.meta.get('SlideXPosUm')

    @property
    def SlideYPosUm(self) -> Optional[str]:
        return self.meta.get('SlideYPosUm')

    @property
    def PanoramaPixelXPos(self) -> Optional[str]:
        return self.meta.get('PanoramaPixelXPos')

    @property
    def PanoramaPixelYPos(self) -> Optional[str]:
        return self.meta.get('PanoramaPixelYPos')

    def __repr__(self):
        return f"<ROIPoint(id={self.id}, roi_id={self.roi_id}, metaname={self.metaname})>"
