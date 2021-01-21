import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB

from histocat.core.base import Base
from histocat.core.utils import remove_location_upon_delete


@remove_location_upon_delete
class SlideModel(Base):
    """Slide."""

    __tablename__ = "slide"
    __table_args__ = (sa.UniqueConstraint("project_id", "name", name="uq_slide_project_id_and_name"),)

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    project_id = sa.Column(sa.Integer(), sa.ForeignKey("project.id", ondelete="CASCADE"), index=True, nullable=False)
    origin_id = sa.Column(sa.Integer(), index=True)
    name = sa.Column(sa.String(4096))
    width_um = sa.Column(sa.Integer())
    height_um = sa.Column(sa.Integer())
    has_slide_image = sa.Column(sa.Boolean(), nullable=False, default=False)
    meta = sa.Column(JSONB())
    session_meta = sa.Column(JSONB())
    location = sa.Column(sa.String(4096))
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    project = sa.orm.relationship("ProjectModel", back_populates="slides")
    panoramas = sa.orm.relationship("PanoramaModel", back_populates="slide", cascade="all, delete, delete-orphan")
    acquisitions = sa.orm.relationship("AcquisitionModel", back_populates="slide", cascade="all, delete, delete-orphan")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
