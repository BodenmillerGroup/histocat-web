import os
import logging
from datetime import datetime

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship

from app.core.utils import remove_location_upon_delete, autocreate_directory_property
from app.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for panorama locations
PANORAMA_LOCATION_FORMAT = "panorama_{id}"


@remove_location_upon_delete
class Panorama(Base):
    """
    Panorama
    """

    __tablename__ = "panorama"

    id: int = sa.Column(sa.Integer(), primary_key=True, index=True)
    slide_id: int = sa.Column(sa.Integer(), sa.ForeignKey("slide.id", ondelete="CASCADE"), index=True)
    metaname: str = sa.Column('metaname', sa.String(4096))
    original_id: int = sa.Column('original_id', sa.Integer(), index=True)
    description: str = sa.Column('description', sa.String())
    slide_x1_pos_um: float = sa.Column('slide_x1_pos_um', sa.Float())
    slide_y1_pos_um: float = sa.Column('slide_y1_pos_um', sa.Float())
    slide_x2_pos_um: float = sa.Column('slide_x2_pos_um', sa.Float())
    slide_y2_pos_um: float = sa.Column('slide_y2_pos_um', sa.Float())
    slide_x3_pos_um: float = sa.Column('slide_x3_pos_um', sa.Float())
    slide_y3_pos_um: float = sa.Column('slide_y3_pos_um', sa.Float())
    slide_x4_pos_um: float = sa.Column('slide_x4_pos_um', sa.Float())
    slide_y4_pos_um: float = sa.Column('slide_y4_pos_um', sa.Float())
    image_end_offset: int = sa.Column('image_end_offset', sa.BigInteger())
    image_start_offset: int = sa.Column('image_start_offset', sa.BigInteger())
    pixel_width: int = sa.Column('pixel_width', sa.Integer())
    pixel_height: int = sa.Column('pixel_height', sa.Integer())
    image_format: str = sa.Column('image_format', sa.String(16))
    pixel_scale_coef: float = sa.Column('pixel_scale_coef', sa.Float())
    meta: dict = sa.Column('meta', JSONB())
    location: str = sa.Column('location', sa.String(4096))
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    slide = relationship("Slide", back_populates="panoramas")
    rois = relationship("ROI", back_populates="panorama", cascade="all, delete, delete-orphan")

    @autocreate_directory_property
    def rois_location(self) -> str:
        """
        Location where ROIs are stored
        """
        return os.path.join(self.location, "rois")

    def __repr__(self):
        return f"<Panorama(id={self.id}, description={self.description})>"
