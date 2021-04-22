import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from histocat.core.base import Base


class GateModel(Base):
    """Cell gate."""

    __tablename__ = "gate"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    dataset_id = sa.Column(sa.Integer(), sa.ForeignKey("dataset.id", ondelete="CASCADE"), index=True, nullable=False)
    name = sa.Column(sa.String())
    description = sa.Column(sa.String())
    cell_classes = sa.Column("cell_classes", JSONB())
    annotations = sa.Column("annotations", JSONB())
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    dataset = sa.orm.relationship("DatasetModel", back_populates="gates")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
