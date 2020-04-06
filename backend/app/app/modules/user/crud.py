from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from app.core.security import get_password_hash, verify_password

from .db import User
from .models import UserCreateModel, UserUpdateModel


def get(session: Session, *, id: int) -> Optional[User]:
    return session.query(User).filter(User.id == id).first()


def get_by_email(session: Session, *, email: str) -> Optional[User]:
    return session.query(User).filter(User.email == email).first()


def authenticate(session: Session, *, email: str, password: str) -> Optional[User]:
    user = get_by_email(session, email=email)
    if not user:
        return None
    if not verify_password(password, user.password):
        return None
    return user


def get_multi(session: Session, *, skip=0, limit=1000) -> List[Optional[User]]:
    return session.query(User).offset(skip).limit(limit).all()


def create(session: Session, *, params: UserCreateModel) -> User:
    entity = User(
        email=params.email,
        password=get_password_hash(params.password),
        name=params.name,
        is_admin=params.is_admin,
        is_active=params.is_active,
    )
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def update(session: Session, *, item: User, params: UserUpdateModel) -> User:
    data = jsonable_encoder(item)
    for field in data:
        if field in params.fields:
            value_in = getattr(params, field)
            if value_in is not None:
                setattr(item, field, value_in)
    if params.password:
        passwordhash = get_password_hash(params.password)
        item.password = passwordhash
    session.add(item)
    session.commit()
    session.refresh(item)
    return item
