import logging
from datetime import datetime
from typing import Optional

import sqlalchemy as sa
from imctools.io import mcdxmlparser
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship

from app.db.base import Base

logger = logging.getLogger(__name__)


class ROI(Base):
    """
    ROI
    """

    __tablename__ = "roi"

    id: int = sa.Column(sa.Integer(), primary_key=True, index=True)
    panorama_id: int = sa.Column(
        sa.Integer(), sa.ForeignKey("panorama.id", ondelete="CASCADE"), index=True
    )
    origin_id: int = sa.Column("origin_id", sa.Integer(), index=True)
    meta: dict = sa.Column("meta", JSONB())
    created_at: datetime = sa.Column(
        "created_at", sa.DateTime(), default=sa.sql.func.now(), nullable=False
    )

    panorama = relationship("Panorama", back_populates="rois")
    acquisitions = relationship(
        "Acquisition", back_populates="roi", cascade="all, delete, delete-orphan"
    )
    roi_points = relationship(
        "ROIPoint", back_populates="roi", cascade="all, delete, delete-orphan"
    )

    @property
    def ROIType(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ROITYPE)

    def __repr__(self):
        return f"<ROI(id={self.id})>"
