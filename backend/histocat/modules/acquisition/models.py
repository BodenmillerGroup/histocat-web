import logging

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from histocat.db.base import Base

logger = logging.getLogger(__name__)


class AcquisitionModel(Base):
    """
    An *acquisition* contains all files belonging to one microscope image acquisition process.
    """

    __tablename__ = "acquisition"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    slide_id = sa.Column(sa.Integer(), sa.ForeignKey("slide.id", ondelete="CASCADE"), index=True, nullable=False)
    origin_id = sa.Column(sa.Integer(), index=True)
    description = sa.Column(sa.String())
    max_x = sa.Column(sa.Integer())
    max_y = sa.Column(sa.Integer())
    signal_type = sa.Column(sa.String())
    segment_data_format = sa.Column(sa.String())
    ablation_frequency = sa.Column(sa.Float())
    ablation_power = sa.Column(sa.Float())
    start_timestamp = sa.Column(sa.DateTime(timezone=True))
    end_timestamp = sa.Column(sa.DateTime(timezone=True))
    movement_type = sa.Column(sa.String())
    ablation_distance_between_shots_x = sa.Column(sa.Float())
    ablation_distance_between_shots_y = sa.Column(sa.Float())
    template = sa.Column(sa.String())
    roi_start_x_pos_um = sa.Column(sa.Float())
    roi_start_y_pos_um = sa.Column(sa.Float())
    roi_end_x_pos_um = sa.Column(sa.Float())
    roi_end_y_pos_um = sa.Column(sa.Float())
    has_before_ablation_image = sa.Column(sa.Boolean(), nullable=False, default=False)
    has_after_ablation_image = sa.Column(sa.Boolean(), nullable=False, default=False)
    is_valid = sa.Column(sa.Boolean(), nullable=False, default=True)
    channels = sa.Column(JSONB())
    meta = sa.Column(JSONB())
    location = sa.Column(sa.String(4096))

    slide = sa.orm.relationship("SlideModel", back_populates="acquisitions")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id})>"
