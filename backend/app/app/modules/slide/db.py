import logging
import os
from datetime import datetime
from typing import Dict, Optional

import sqlalchemy as sa
from imctools.io import mcdxmlparser
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
    __table_args__ = (
        sa.UniqueConstraint('experiment_id', 'uid', name='uq_experiment_slide_uid'),
    )

    id: int = sa.Column(sa.Integer(), primary_key=True, index=True)
    experiment_id: int = sa.Column(sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), index=True)
    metaname: str = sa.Column('metaname', sa.String(4096))
    original_id: int = sa.Column('original_id', sa.Integer(), index=True)
    uid: str = sa.Column('uid', UUID(), index=True)
    original_metadata: str = sa.Column('original_metadata', sa.String())
    meta: Dict[str, str] = sa.Column('meta', JSONB())
    location: str = sa.Column('location', sa.String(4096))
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    experiment = relationship("Experiment", back_populates="slides")
    panoramas = relationship("Panorama", back_populates="slide", cascade="all, delete, delete-orphan")

    @property
    def UID(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.UID)

    @property
    def Description(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.DESCRIPTION)

    @property
    def Filename(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.FILENAME)

    @property
    def SlideType(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.SLIDETYPE)

    @property
    def WidthUm(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.WIDTHUM)

    @property
    def HeightUm(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.HEIGHTUM)

    @property
    def ImageEndOffset(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.IMAGEENDOFFSET)

    @property
    def ImageStartOffset(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.IMAGESTARTOFFSET)

    @property
    def ImageFile(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.IMAGEFILE)

    @autocreate_directory_property
    def panoramas_location(self) -> str:
        """
        Location where panoramas are stored
        """
        return os.path.join(self.location, "panoramas")

    def __repr__(self):
        return f"<Slide(id={self.id}, metaname={self.metaname})>"
