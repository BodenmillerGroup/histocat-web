from typing import List, Optional

from fastapi.encoders import jsonable_encoder

from app.core.security import get_password_hash, verify_password
from app.modules.user.db import User
from .models import UserInCreateModel, UserInUpdateModel


def get(db_session, *, id: int) -> Optional[User]:
    return db_session.query(User).filter(User.id == id).first()


def get_by_email(db_session, *, email: str) -> Optional[User]:
    return db_session.query(User).filter(User.email == email).first()


def authenticate(db_session, *, email: str, password: str) -> Optional[User]:
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


def get_multi(db_session, *, skip=0, limit=100) -> List[Optional[User]]:
    return db_session.query(User).offset(skip).limit(limit).all()


def create(db_session, *, params: UserInCreateModel) -> User:
    user = User(
        email=params.email,
        hashed_password=get_password_hash(params.password),
        full_name=params.full_name,
        is_superuser=params.is_superuser,
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    return user


def update(db_session, *, user: User, params: UserInUpdateModel) -> User:
    user_data = jsonable_encoder(user)
    for field in user_data:
        if field in params.fields:
            value_in = getattr(params, field)
            if value_in is not None:
                setattr(user, field, value_in)
    if params.password:
        passwordhash = get_password_hash(params.password)
        user.hashed_password = passwordhash
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    return user
