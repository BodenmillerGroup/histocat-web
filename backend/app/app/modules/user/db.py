from datetime import datetime

from sqlalchemy import Boolean, Column, String, Integer, DateTime
from sqlalchemy.orm import relationship
from sqlalchemy.sql.functions import now

from app.db.base import Base


class User(Base):
    __tablename__ = "user"

    id: int = Column(Integer, primary_key=True, index=True)
    full_name: str = Column(String, index=True)
    email: str = Column(String, unique=True, index=True, nullable=False)
    hashed_password: str = Column(String)
    is_active: bool = Column(Boolean, default=True, nullable=False)
    is_superuser: bool = Column(Boolean, default=False, nullable=False)
    created_at: datetime = Column(DateTime, default=now(), nullable=False)
    updated_at: datetime = Column(DateTime, default=now(), onupdate=now(), nullable=False)

    experiments = relationship("Experiment", back_populates="user", cascade="all, delete, delete-orphan")

    def __repr__(self):
        return f"<User(id={self.id}, full_name={self.full_name}, email={self.email}, is_active={self.is_active}, is_superuser={self.is_superuser})>"
