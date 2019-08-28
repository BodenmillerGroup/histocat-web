import logging
from datetime import datetime

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship

from app.db.base import Base

logger = logging.getLogger(__name__)


class AcquisitionArtifact(Base):
    """
    Acquisition artifact (mask, image, csv file, etc.
    """

    __tablename__ = "acquisition_artifact"

    id: int = sa.Column('id', sa.Integer(), primary_key=True, index=True)
    acquisition_id: int = sa.Column('acquisition_id', sa.Integer(), sa.ForeignKey("acquisition.id", ondelete="CASCADE"),
                                    index=True)
    type: str = sa.Column('type', sa.String(64), nullable=False, index=True)
    description: str = sa.Column('description', sa.String())
    meta: dict = sa.Column('meta', JSONB())
    location: str = sa.Column('location', sa.String(4096))
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    acquisition = relationship("Acquisition", back_populates="artifacts")

    def __repr__(self):
        return f"<AcquisitionArtifact(id={self.id}, type={self.type}, description={self.description})>"
