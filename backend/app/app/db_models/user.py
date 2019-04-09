from datetime import datetime

from sqlalchemy import Boolean, Column, Integer, String, DateTime

from app.db.base_class import Base


class User(Base):
    id = Column(Integer, primary_key=True, index=True)
    full_name = Column(String, index=True)
    email = Column(String, unique=True, index=True)
    hashed_password = Column(String)
    is_active = Column(Boolean, default=True)
    is_superuser = Column(Boolean, default=False)
    created_at = Column(DateTime, default=datetime.utcnow)

    def __repr__(self):
        return f'<User(id={self.id}, full_name={self.full_name}, email={self.email}, is_active={self.is_active},' \
            f' is_superuser={self.is_superuser})>'
