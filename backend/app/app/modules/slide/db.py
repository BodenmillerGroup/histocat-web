import logging
import os
from datetime import datetime

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB, UUID
from sqlalchemy.orm import relationship

from app.core.utils import autocreate_directory_property, remove_location_upon_delete
from app.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for slide locations
SLIDE_LOCATION_FORMAT = "slide_{id}"


@remove_location_upon_delete
class Slide(Base):
    """
    Slide
    """

    __tablename__ = "slide"

    id: int = sa.Column(sa.Integer(), primary_key=True, index=True)
    experiment_id: int = sa.Column(sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), index=True)
    uid: str = sa.Column('uid', UUID(), index=True)
    description: str = sa.Column('description', sa.String())
    filename: str = sa.Column('filename', sa.String(4096))
    slide_type: str = sa.Column('slide_type', sa.String())
    width_um: int = sa.Column('width_um', sa.Integer())
    height_um: int = sa.Column('height_um', sa.Integer())
    image_end_offset: int = sa.Column('image_end_offset', sa.BigInteger())
    image_start_offset: int = sa.Column('image_start_offset', sa.BigInteger())
    image_file: str = sa.Column('image_file', sa.String(4096))
    meta: dict = sa.Column('meta', JSONB())
    location: str = sa.Column('location', sa.String(4096))
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    experiment = relationship("Experiment", back_populates="slides")
    panoramas = relationship("Panorama", back_populates="slide", cascade="all, delete, delete-orphan")

    @autocreate_directory_property
    def panoramas_location(self) -> str:
        """
        Location where panoramas are stored
        """
        return os.path.join(self.location, "panoramas")

    def __repr__(self):
        return f"<Slide(id={self.id}, description={self.description})>"
