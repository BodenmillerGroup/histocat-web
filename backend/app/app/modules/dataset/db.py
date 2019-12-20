import logging
from datetime import datetime

import sqlalchemy as sa
from sqlalchemy import text
from sqlalchemy.dialects.postgresql import JSONB, UUID
from sqlalchemy.orm import relationship

from app.core.utils import remove_location_upon_delete
from app.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for dataset locations
DATASET_LOCATION_FORMAT = "dataset_{id}"


@remove_location_upon_delete
class Dataset(Base):
    """
    Dataset
    """

    __tablename__ = "dataset"

    id: int = sa.Column("id", sa.Integer(), primary_key=True, index=True)
    experiment_id: int = sa.Column(
        "experiment_id", sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), index=True,
    )
    user_id: int = sa.Column(
        "user_id", sa.Integer(), sa.ForeignKey("user.id", ondelete="CASCADE"), index=True,
    )
    uid: str = sa.Column(
        "uid", UUID(), server_default=text("uuid_generate_v4()"), nullable=False, index=True,
    )
    status: str = sa.Column("status", sa.String(64), default="pending", nullable=False, index=True)
    name: str = sa.Column("name", sa.String())
    description: str = sa.Column("description", sa.String())
    input: dict = sa.Column("input", JSONB())
    output: dict = sa.Column("output", JSONB())
    meta: dict = sa.Column("meta", JSONB())
    location: str = sa.Column("location", sa.String(4096))
    created_at: datetime = sa.Column("created_at", sa.DateTime(), default=sa.sql.func.now(), nullable=False)
    updated_at: datetime = sa.Column(
        sa.DateTime(), default=sa.sql.func.now(), onupdate=sa.sql.func.now(), nullable=False,
    )

    experiment = relationship("Experiment", back_populates="datasets")
    user = relationship("User", back_populates="datasets")

    def __repr__(self):
        return f"<Dataset(id={self.id}, name={self.name})>"
