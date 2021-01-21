import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from histocat.core.base import Base


class PipelineModel(Base):
    """Project processing pipeline."""

    __tablename__ = "pipeline"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    project_id = sa.Column(sa.Integer(), sa.ForeignKey("project.id", ondelete="CASCADE"), index=True, nullable=False)
    name = sa.Column(sa.String())
    description = sa.Column(sa.String())
    steps = sa.Column(JSONB())
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    project = sa.orm.relationship("ProjectModel", back_populates="pipelines")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
