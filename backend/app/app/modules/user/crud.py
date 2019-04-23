from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from app.core.security import get_password_hash, verify_password
from .db import User
from .models import UserCreateModel, UserUpdateModel


def get(db_session: Session, *, id: int) -> Optional[User]:
    return db_session.query(User).filter(User.id == id).first()


def get_by_email(db_session: Session, *, email: str) -> Optional[User]:
    return db_session.query(User).filter(User.email == email).first()


def authenticate(db_session: Session, *, email: str, password: str) -> Optional[User]:
    user = get_by_email(db_session, email=email)
    if not user:
        return None
    if not verify_password(password, user.hashed_password):
        return None
    return user


def is_active(user) -> bool:
    return user.is_active


def is_superuser(user) -> bool:
    return user.is_superuser


def get_multi(db_session: Session, *, skip=0, limit=100) -> List[Optional[User]]:
    return db_session.query(User).offset(skip).limit(limit).all()


def create(db_session: Session, *, params: UserCreateModel) -> User:
    entity = User(
        email=params.email,
        hashed_password=get_password_hash(params.password),
        full_name=params.full_name,
        is_superuser=params.is_superuser,
    )
    db_session.add(entity)
    db_session.commit()
    db_session.refresh(entity)
    return entity


def update(db_session: Session, *, item: User, params: UserUpdateModel) -> User:
    data = jsonable_encoder(item)
    for field in data:
        if field in params.fields:
            value_in = getattr(params, field)
            if value_in is not None:
                setattr(item, field, value_in)
    if params.password:
        passwordhash = get_password_hash(params.password)
        item.hashed_password = passwordhash
    db_session.add(item)
    db_session.commit()
    db_session.refresh(item)
    return item
