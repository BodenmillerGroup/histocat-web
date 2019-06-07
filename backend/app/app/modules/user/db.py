from sqlalchemy import Boolean, Column, String

from app.db.base import Base, CreatedAtMixin, IdMixin, UpdatedAtMixin


class User(IdMixin, Base, CreatedAtMixin, UpdatedAtMixin):
    __tablename__ = "user"

    full_name = Column(String, index=True)
    email = Column(String, unique=True, index=True, nullable=False)
    hashed_password = Column(String)
    is_active = Column(Boolean, default=True, nullable=False)
    is_superuser = Column(Boolean, default=False, nullable=False)

    def __repr__(self):
        return (
            f"<User(id={self.id}, full_name={self.full_name}, email={self.email}, is_active={self.is_active}, is_superuser={self.is_superuser})>"
        )
