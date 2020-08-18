import sqlalchemy as sa

from histocat.db.base import Base


class GroupModel(Base):
    __tablename__ = "group"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    name = sa.Column(sa.String(), index=True, unique=True, nullable=False)
    description = sa.Column(sa.String())
    url = sa.Column(sa.String())
    is_open = sa.Column(sa.Boolean(), default=False, nullable=False)
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)

    members = sa.orm.relationship("MemberModel", back_populates="group", cascade="all, delete, delete-orphan")
    experiments = sa.orm.relationship("ExperimentModel", back_populates="group", cascade="all, delete, delete-orphan")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, name={self.name})>"
