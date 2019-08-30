import logging
from datetime import datetime
from typing import Optional

import sqlalchemy as sa
from imctools.io import mcdxmlparser
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
    origin_id: int = sa.Column('origin_id', sa.Integer(), index=True)
    meta: dict = sa.Column('meta', JSONB())
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    roi = relationship("ROI", back_populates="roi_points")

    @property
    def OrderNumber(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ORDERNUMBER)

    @property
    def SlideXPosUm(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.SLIDEXPOSUM)

    @property
    def SlideYPosUm(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.SLIDEYPOSUM)

    @property
    def PanoramaPixelXPos(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.PANORAMAPIXELXPOS)

    @property
    def PanoramaPixelYPos(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.PANORAMAPIXELYPOS)

    def __repr__(self):
        return f"<ROIPoint(id={self.id})>"
