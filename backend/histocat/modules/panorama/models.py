import logging

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from histocat.db.base import Base

logger = logging.getLogger(__name__)


class PanoramaModel(Base):
    """Panoramic image."""

    __tablename__ = "panorama"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    slide_id = sa.Column(sa.Integer(), sa.ForeignKey("slide.id", ondelete="CASCADE"), index=True)
    origin_id = sa.Column(sa.Integer(), index=True)
    image_type = sa.Column(sa.String())
    description = sa.Column(sa.String())
    start_position_x = sa.Column(sa.Float())
    start_position_y = sa.Column(sa.Float())
    width = sa.Column(sa.Float())
    height = sa.Column(sa.Float())
    rotation_angle = sa.Column(sa.Float())
    meta = sa.Column(JSONB())
    location = sa.Column(sa.String(4096))

    slide = sa.orm.relationship("SlideModel", back_populates="panoramas")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id})>"
