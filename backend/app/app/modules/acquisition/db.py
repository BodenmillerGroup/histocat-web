import logging
import os
from datetime import datetime

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
    description: str = sa.Column('description', sa.String())
    metaname: str = sa.Column('metaname', sa.String(4096))
    original_id: int = sa.Column('original_id', sa.Integer(), index=True)
    order_number: int = sa.Column('order_number', sa.Integer())
    ablation_power: float = sa.Column('ablation_power', sa.Float())
    ablation_distance_between_shots_x: float = sa.Column('ablation_distance_between_shots_x', sa.Float())
    ablation_distance_between_shots_y: float= sa.Column('ablation_distance_between_shots_y', sa.Float())
    ablation_frequency: float = sa.Column('ablation_frequency', sa.Float())
    signal_type: str = sa.Column('signal_type', sa.String(64))
    dual_count_start: int = sa.Column('dual_count_start', sa.Integer())
    data_start_offset: int = sa.Column('data_start_offset', sa.BigInteger())
    data_end_offset: int = sa.Column('data_end_offset', sa.BigInteger())
    start_timestamp: datetime = sa.Column('start_timestamp', sa.DateTime())
    end_timestamp: datetime = sa.Column('end_timestamp', sa.DateTime())
    after_ablation_image_end_offset: int = sa.Column('after_ablation_image_end_offset', sa.BigInteger())
    after_ablation_image_start_offset: int = sa.Column('after_ablation_image_start_offset', sa.BigInteger())
    before_ablation_image_end_offset: int = sa.Column('before_ablation_image_end_offset', sa.BigInteger())
    before_ablation_image_start_offset: int = sa.Column('before_ablation_image_start_offset', sa.BigInteger())
    roi_start_x_pos_um: float = sa.Column('roi_start_x_pos_um', sa.Float())
    roi_start_y_pos_um: float = sa.Column('roi_start_y_pos_um', sa.Float())
    roi_end_x_pos_um: float = sa.Column('roi_end_x_pos_um', sa.Float())
    roi_end_y_pos_um: float = sa.Column('roi_end_y_pos_um', sa.Float())
    movement_type: str = sa.Column('movement_type', sa.String(64))
    segment_data_format: str = sa.Column('segment_data_format', sa.String(64))
    value_bytes: int = sa.Column('value_bytes', sa.Integer())
    max_y: int = sa.Column('max_y', sa.Integer())
    max_x: int = sa.Column('max_x', sa.Integer())
    plume_start: int = sa.Column('plume_start', sa.Integer())
    plume_end: int = sa.Column('plume_end', sa.Integer())
    template: str = sa.Column('template', sa.String())
    meta: dict = sa.Column('meta', JSONB())
    location: str = sa.Column('location', sa.String(4096))
    created_at: datetime = sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    roi = relationship("ROI", back_populates="acquisitions")
    channels = relationship("Channel", back_populates="acquisition", cascade="all, delete, delete-orphan")

    @autocreate_directory_property
    def channels_location(self) -> str:
        """str: location where channels files are stored"""
        return os.path.join(self.location, "channels")

    def __repr__(self):
        return f"<Acquisition(id={self.id}, description={self.description})>"
