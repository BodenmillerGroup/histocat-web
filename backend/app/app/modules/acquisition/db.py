import logging
from datetime import datetime
from typing import Optional

import sqlalchemy as sa
from imctools.io import mcdxmlparser
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship

from app.db.base import Base

logger = logging.getLogger(__name__)


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
        return self.meta.get(mcdxmlparser.DESCRIPTION)

    @property
    def OrderNumber(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ORDERNUMBER)

    @property
    def AblationPower(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ABLATIONPOWER)

    @property
    def AblationDistanceBetweenShotsX(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ABLATIONDISTANCEBETWEENSHOTSX)

    @property
    def AblationDistanceBetweenShotsY(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ABLATIONDISTANCEBETWEENSHOTSY)

    @property
    def AblationFrequency(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ABLATIONFREQUENCY)

    @property
    def SignalType(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.SIGNALTYPE)

    @property
    def DualCountStart(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.DUALCOUNTSTART)

    @property
    def DataStartOffset(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.DATASTARTOFFSET)

    @property
    def DataEndOffset(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.DATAENDOFFSET)

    @property
    def StartTimeStamp(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.STARTTIMESTAMP)

    @property
    def EndTimeStamp(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ENDTIMESTAMP)

    @property
    def AfterAblationImageEndOffset(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.AFTERABLATIONIMAGEENDOFFSET)

    @property
    def AfterAblationImageStartOffset(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.AFTERABLATIONIMAGESTARTOFFSET)

    @property
    def BeforeAblationImageEndOffset(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.BEFOREABLATIONIMAGEENDOFFSET)

    @property
    def BeforeAblationImageStartOffset(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.BEFOREABLATIONIMAGESTARTOFFSET)

    @property
    def ROIStartXPosUm(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ROISTARTXPOSUM)

    @property
    def ROIStartYPosUm(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ROISTARTYPOSUM)

    @property
    def ROIEndXPosUm(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ROIENDXPOSUM)

    @property
    def ROIEndYPosUm(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.ROIENDYPOSUM)

    @property
    def MovementType(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.MOVEMENTTYPE)

    @property
    def SegmentDataFormat(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.SEGMENTDATAFORMAT)

    @property
    def ValueBytes(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.VALUEBYTES)

    @property
    def MaxY(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.MAXY)

    @property
    def MaxX(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.MAXX)

    @property
    def PlumeStart(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.PLUMESTART)

    @property
    def PlumeEnd(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.PLUMEEND)

    @property
    def Template(self) -> Optional[str]:
        return self.meta.get(mcdxmlparser.TEMPLATE)

    def __repr__(self):
        return f"<Acquisition(id={self.id}, metaname={self.metaname})>"
