import os
import logging
from datetime import datetime

import sqlalchemy as sa
from sqlalchemy.orm import relationship

from app.core.utils import remove_location_upon_delete, autocreate_directory_property
from app.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for ROI locations
ROI_LOCATION_FORMAT = "roi_{id}"


@remove_location_upon_delete
class ROI(Base):
    """
    ROI
    """

    __tablename__ = "roi"

    id: int = sa.Column(sa.Integer(), primary_key=True, index=True)
    panorama_id: int = sa.Column(sa.Integer(), sa.ForeignKey("panorama.id", ondelete="CASCADE"), index=True)
    roi_type: str = sa.Column('roi_type', sa.String(64))
    location: str = sa.Column('location', sa.String(4096))
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    panorama = relationship("Panorama", back_populates="rois")
    acquisitions = relationship("Acquisition", back_populates="roi", cascade="all, delete, delete-orphan")
    roi_points = relationship("ROIPoint", back_populates="roi", cascade="all, delete, delete-orphan")

    @autocreate_directory_property
    def acquisitions_location(self) -> str:
        """
        Location where acquisitions are stored
        """
        return os.path.join(self.location, "acquisitions")

    def __repr__(self):
        return f"<ROI(id={self.id}, roi_type={self.roi_type})>"
