from sqlalchemy import Boolean, Column, Integer, String, DateTime, UniqueConstraint
from sqlalchemy.sql.functions import now

from app.db_models.base import Base


class User(Base):
    __tablename__ = 'user'

    id = Column(Integer, primary_key=True, autoincrement=True, index=True)
    full_name = Column(String, index=True)
    email = Column(String, unique=True, index=True, nullable=False)
    hashed_password = Column(String)
    is_active = Column(Boolean, default=True)
    is_superuser = Column(Boolean, default=False)
    created_at = Column(DateTime, default=now())

    def __repr__(self):
        return f'<User(id={self.id}, full_name={self.full_name}, email={self.email}, is_active={self.is_active},' \
            f' is_superuser={self.is_superuser})>'
