from datetime import datetime

import sqlalchemy as sa
from sqlalchemy.orm import relationship

from app.db.base import Base


class User(Base):
    __tablename__ = "user"

    id: int = sa.Column(sa.Integer(), primary_key=True, index=True)
    full_name: str = sa.Column(sa.String(), index=True)
    email: str = sa.Column(sa.String(), unique=True, index=True, nullable=False)
    hashed_password: str = sa.Column(sa.String())
    is_active: bool = sa.Column(sa.Boolean(), default=True, nullable=False)
    is_superuser: bool = sa.Column(sa.Boolean(), default=False, nullable=False)
    created_at: datetime = sa.Column(sa.DateTime(), default=sa.sql.func.now(), nullable=False)
    updated_at: datetime = sa.Column(
        sa.DateTime(), default=sa.sql.func.now(), onupdate=sa.sql.func.now(), nullable=False,
    )

    experiments = relationship("Experiment", back_populates="user", cascade="all, delete, delete-orphan")
    datasets = relationship("Dataset", back_populates="user", cascade="all, delete, delete-orphan")

    def __repr__(self):
        return f"<User(id={self.id}, full_name={self.full_name}, email={self.email}, is_active={self.is_active}, is_superuser={self.is_superuser})>"
