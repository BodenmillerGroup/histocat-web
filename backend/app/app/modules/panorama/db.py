import logging
import os
from datetime import datetime
from typing import Optional

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
    meta: dict = sa.Column('meta', JSONB())
    location: str = sa.Column('location', sa.String(4096))
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    slide = relationship("Slide", back_populates="panoramas")
    rois = relationship("ROI", back_populates="panorama", cascade="all, delete, delete-orphan")

    @property
    def Description(self) -> Optional[str]:
        return self.meta.get('Description')

    @property
    def SlideX1PosUm(self) -> Optional[str]:
        return self.meta.get('SlideX1PosUm')

    @property
    def SlideY1PosUm(self) -> Optional[str]:
        return self.meta.get('SlideY1PosUm')

    @property
    def SlideX2PosUm(self) -> Optional[str]:
        return self.meta.get('SlideX2PosUm')

    @property
    def SlideY2PosUm(self) -> Optional[str]:
        return self.meta.get('SlideY2PosUm')

    @property
    def SlideX3PosUm(self) -> Optional[str]:
        return self.meta.get('SlideX3PosUm')

    @property
    def SlideY3PosUm(self) -> Optional[str]:
        return self.meta.get('SlideY3PosUm')

    @property
    def SlideX4PosUm(self) -> Optional[str]:
        return self.meta.get('SlideX4PosUm')

    @property
    def SlideY4PosUm(self) -> Optional[str]:
        return self.meta.get('SlideY4PosUm')

    @property
    def ImageEndOffset(self) -> Optional[str]:
        return self.meta.get('ImageEndOffset')

    @property
    def ImageStartOffset(self) -> Optional[str]:
        return self.meta.get('ImageStartOffset')

    @property
    def PixelWidth(self) -> Optional[str]:
        return self.meta.get('PixelWidth')

    @property
    def PixelHeight(self) -> Optional[str]:
        return self.meta.get('PixelHeight')

    @property
    def ImageFormat(self) -> Optional[str]:
        return self.meta.get('ImageFormat')

    @property
    def PixelScaleCoef(self) -> Optional[str]:
        return self.meta.get('PixelScaleCoef')

    @autocreate_directory_property
    def rois_location(self) -> str:
        """
        Location where ROIs are stored
        """
        return os.path.join(self.location, "rois")

    def __repr__(self):
        return f"<Panorama(id={self.id}, metaname={self.metaname})>"
