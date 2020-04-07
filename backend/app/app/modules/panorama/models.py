import logging

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from app.db.base import Base

logger = logging.getLogger(__name__)


class Panorama(Base):
    """Panoramic image."""

    __tablename__ = "panorama"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    slide_id = sa.Column(sa.Integer(), sa.ForeignKey("slide.id", ondelete="CASCADE"), index=True)
    origin_id = sa.Column(sa.Integer(), index=True)
    meta = sa.Column(JSONB())
    location = sa.Column(sa.String(4096))
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    slide = sa.orm.relationship("Slide", back_populates="panoramas")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id})>"
