import logging

import sqlalchemy as sa
from sqlalchemy import text
from sqlalchemy.dialects.postgresql import JSONB, UUID

from histocat.core.utils import remove_location_upon_delete
from histocat.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for dataset locations
DATASET_LOCATION_FORMAT = "dataset_{id}"


@remove_location_upon_delete
class DatasetModel(Base):
    """Dataset."""

    __tablename__ = "dataset"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    experiment_id = sa.Column(sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), index=True)
    user_id = sa.Column(sa.Integer(), sa.ForeignKey("user.id", ondelete="CASCADE"), index=True)
    uid = sa.Column(UUID(), server_default=text("uuid_generate_v4()"), nullable=False, index=True)
    status = sa.Column(sa.String(64), default="pending", nullable=False, index=True)
    name = sa.Column(sa.String())
    description = sa.Column(sa.String())
    input = sa.Column(JSONB())
    output = sa.Column(JSONB())
    meta = sa.Column(JSONB())
    location = sa.Column(sa.String(4096))
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)
    updated_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), onupdate=sa.sql.func.now(), nullable=False)

    experiment = sa.orm.relationship("ExperimentModel", back_populates="datasets")
    user = sa.orm.relationship("UserModel", back_populates="datasets")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
