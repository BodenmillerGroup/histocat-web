from typing import Optional, Sequence

from histocat.core.security import verify_password
from sqlalchemy.orm import Session

from .dto import UserCreateDto, UserUpdateDto
from .models import UserModel


def get_by_id(session: Session, id: int) -> Optional[UserModel]:
    return session.query(UserModel).filter(UserModel.id == id).first()


def get_by_email(session: Session, *, email: str) -> Optional[UserModel]:
    return session.query(UserModel).filter(UserModel.email == email).first()


def authenticate(session: Session, *, email: str, password: str) -> Optional[UserModel]:
    user = get_by_email(session, email=email)
    if not user:
        return None
    if not verify_password(password, user.password):
        return None
    return user


def get_all(session: Session, skip=0, limit=1000) -> Sequence[UserModel]:
    return session.query(UserModel).offset(skip).limit(limit).all()


def create(session: Session, *, params: UserCreateDto) -> UserModel:
    entity = UserModel(
        email=params.email,
        password=params.password,
        name=params.name,
        is_admin=params.is_admin,
        is_active=params.is_active,
    )
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def update(session: Session, *, item: UserModel, params: UserUpdateDto) -> UserModel:
    data = item.as_dict()
    update_data = params.dict(exclude_unset=True)
    for field in data:
        if field in update_data:
            setattr(item, field, update_data[field])
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def confirm_signup(session: Session, *, item: UserModel):
    item.is_active = True
    session.add(item)
    session.commit()
    session.refresh(item)
    return item
