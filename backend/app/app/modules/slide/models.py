import logging

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from app.core.utils import remove_location_upon_delete
from app.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for slide locations
SLIDE_LOCATION_FORMAT = "slide_{id}"


@remove_location_upon_delete
class Slide(Base):
    """Slide."""

    __tablename__ = "slide"
    __table_args__ = (sa.UniqueConstraint("experiment_id", "name", name="uq_experiment_slide_name"),)

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    origin_id = sa.Column(sa.Integer(), index=True)
    experiment_id = sa.Column(sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), index=True)
    name = sa.Column(sa.String(4096))
    meta = sa.Column(JSONB())
    session_meta = sa.Column(JSONB())
    xml_meta = sa.Column(sa.Text())
    location = sa.Column(sa.String(4096))
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    experiment = sa.orm.relationship("Experiment", back_populates="slides")
    panoramas = sa.orm.relationship("Panorama", back_populates="slide", cascade="all, delete, delete-orphan")
    acquisitions = sa.orm.relationship("Acquisition", back_populates="slide", cascade="all, delete, delete-orphan")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
