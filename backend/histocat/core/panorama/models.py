import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from histocat.core.base import Base


class PanoramaModel(Base):
    """Panoramic image."""

    __tablename__ = "panorama"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    slide_id = sa.Column(sa.Integer(), sa.ForeignKey("slide.id", ondelete="CASCADE"), index=True, nullable=False)
    origin_id = sa.Column(sa.Integer(), index=True)
    image_type = sa.Column(sa.String())
    description = sa.Column(sa.String())
    x1 = sa.Column(sa.Float())
    y1 = sa.Column(sa.Float())
    x2 = sa.Column(sa.Float())
    y2 = sa.Column(sa.Float())
    x3 = sa.Column(sa.Float())
    y3 = sa.Column(sa.Float())
    x4 = sa.Column(sa.Float())
    y4 = sa.Column(sa.Float())
    rotation_angle = sa.Column(sa.Float())
    meta = sa.Column(JSONB())
    location = sa.Column(sa.String(4096))

    slide = sa.orm.relationship("SlideModel", back_populates="panoramas")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id})>"
