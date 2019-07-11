import logging
import os
from datetime import datetime
from typing import Optional

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship

from app.core.utils import autocreate_directory_property, remove_location_upon_delete
from app.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for acquisition locations
ACQUISITION_LOCATION_FORMAT = "acquisition_{id}"


@remove_location_upon_delete
class Acquisition(Base):
    """
    An *acquisition* contains all files belonging to one microscope image acquisition process.
    """

    __tablename__ = "acquisition"

    id: int = sa.Column('id', sa.Integer(), primary_key=True, index=True)
    roi_id: int = sa.Column('roi_id', sa.Integer(), sa.ForeignKey("roi.id", ondelete="CASCADE"), index=True)
    metaname: str = sa.Column('metaname', sa.String(4096))
    original_id: int = sa.Column('original_id', sa.Integer(), index=True)
    meta: dict = sa.Column('meta', JSONB())
    location: str = sa.Column('location', sa.String(4096))
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    roi = relationship("ROI", back_populates="acquisitions")
    channels = relationship("Channel", back_populates="acquisition", cascade="all, delete, delete-orphan")

    @property
    def Description(self) -> Optional[str]:
        return self.meta.get('Description')

    @property
    def OrderNumber(self) -> Optional[str]:
        return self.meta.get('OrderNumber')

    @property
    def AblationPower(self) -> Optional[str]:
        return self.meta.get('AblationPower')

    @property
    def AblationDistanceBetweenShotsX(self) -> Optional[str]:
        return self.meta.get('AblationDistanceBetweenShotsX')

    @property
    def AblationDistanceBetweenShotsY(self) -> Optional[str]:
        return self.meta.get('AblationDistanceBetweenShotsY')

    @property
    def AblationFrequency(self) -> Optional[str]:
        return self.meta.get('AblationFrequency')

    @property
    def SignalType(self) -> Optional[str]:
        return self.meta.get('SignalType')

    @property
    def DualCountStart(self) -> Optional[str]:
        return self.meta.get('DualCountStart')

    @property
    def DataStartOffset(self) -> Optional[str]:
        return self.meta.get('DataStartOffset')

    @property
    def DataEndOffset(self) -> Optional[str]:
        return self.meta.get('DataEndOffset')

    @property
    def StartTimeStamp(self) -> Optional[str]:
        return self.meta.get('StartTimeStamp')

    @property
    def EndTimeStamp(self) -> Optional[str]:
        return self.meta.get('EndTimeStamp')

    @property
    def AfterAblationImageEndOffset(self) -> Optional[str]:
        return self.meta.get('AfterAblationImageEndOffset')

    @property
    def AfterAblationImageStartOffset(self) -> Optional[str]:
        return self.meta.get('AfterAblationImageStartOffset')

    @property
    def BeforeAblationImageEndOffset(self) -> Optional[str]:
        return self.meta.get('BeforeAblationImageEndOffset')

    @property
    def BeforeAblationImageStartOffset(self) -> Optional[str]:
        return self.meta.get('BeforeAblationImageStartOffset')

    @property
    def ROIStartXPosUm(self) -> Optional[str]:
        return self.meta.get('ROIStartXPosUm')

    @property
    def ROIStartYPosUm(self) -> Optional[str]:
        return self.meta.get('ROIStartYPosUm')

    @property
    def ROIEndXPosUm(self) -> Optional[str]:
        return self.meta.get('ROIEndXPosUm')

    @property
    def MovementType(self) -> Optional[str]:
        return self.meta.get('MovementType')

    @property
    def SegmentDataFormat(self) -> Optional[str]:
        return self.meta.get('SegmentDataFormat')

    @property
    def ValueBytes(self) -> Optional[str]:
        return self.meta.get('ValueBytes')

    @property
    def MaxY(self) -> Optional[str]:
        return self.meta.get('MaxY')

    @property
    def MaxX(self) -> Optional[str]:
        return self.meta.get('MaxX')

    @property
    def PlumeStart(self) -> Optional[str]:
        return self.meta.get('PlumeStart')

    @property
    def PlumeEnd(self) -> Optional[str]:
        return self.meta.get('PlumeEnd')

    @property
    def Template(self) -> Optional[str]:
        return self.meta.get('Template')

    @autocreate_directory_property
    def channels_location(self) -> str:
        """str: location where channels files are stored"""
        return os.path.join(self.location, "channels")

    def __repr__(self):
        return f"<Acquisition(id={self.id}, metaname={self.metaname})>"
