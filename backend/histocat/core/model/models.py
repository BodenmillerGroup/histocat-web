import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from histocat.core.base import Base
from histocat.core.utils import remove_location_upon_delete


@remove_location_upon_delete
class ModelModel(Base):
    """Segmentation Keras model."""

    __tablename__ = "model"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    name = sa.Column(sa.String(), index=True, unique=True)
    description = sa.Column(sa.Text())
    location = sa.Column(sa.String(4096))
    meta = sa.Column(JSONB())
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
