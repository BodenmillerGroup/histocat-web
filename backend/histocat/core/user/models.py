import sqlalchemy as sa

from histocat.core.base import Base


class UserModel(Base):
    __tablename__ = "user"

    id = sa.Column(sa.Integer(), primary_key=True, index=True)
    name = sa.Column(sa.String(), index=True)
    email = sa.Column(sa.String(), unique=True, index=True, nullable=False)
    password = sa.Column(sa.String())
    is_active = sa.Column(sa.Boolean(), default=True, nullable=False)
    is_admin = sa.Column(sa.Boolean(), default=False, nullable=False)
    created_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)
    updated_at = sa.Column(sa.DateTime(), default=sa.sql.func.now(), onupdate=sa.sql.func.now(), nullable=False)

    members = sa.orm.relationship("MemberModel", back_populates="user", cascade="all, delete, delete-orphan")

    def __repr__(self):
        return f"<{self.__class__.__name__}(id={self.id}, email={self.email}, name={self.name})>"
