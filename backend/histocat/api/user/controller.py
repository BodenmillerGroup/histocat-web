from typing import Sequence

from fastapi import APIRouter, Body, Depends, HTTPException
from fastapi.encoders import jsonable_encoder
from pydantic import EmailStr
from sqlalchemy.orm import Session
from starlette import status

from histocat.api.db import get_db
from histocat.api.security import get_active_user, get_admin
from histocat.config import config
from histocat.core.user import service
from histocat.core.user.dto import UserCreateDto, UserDto, UserUpdateDto
from histocat.core.user.models import UserModel
from histocat.core.utils import send_new_account_email

router = APIRouter()


@router.get("/users", response_model=Sequence[UserDto])
def get_all(db: Session = Depends(get_db), user: UserModel = Depends(get_active_user)):
    """Get all users."""
    items = service.get_all(db)
    return items


@router.post("/users", response_model=UserDto)
def create(
    params: UserCreateDto,
    db: Session = Depends(get_db),
    user: UserModel = Depends(get_admin),
):
    """Create new user."""
    item = service.get_by_email(db, email=params.email)
    if item:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="The user with this username already exists in the system.",
        )
    item = service.create(db, params=params)
    if config.EMAILS_ENABLED and params.email:
        send_new_account_email(email_to=params.email, username=params.name, password=params.password)
    return item


@router.patch("/users/profile", response_model=UserDto)
def update_me(
    password: str = Body(None),
    name: str = Body(None),
    email: EmailStr = Body(None),
    db: Session = Depends(get_db),
    user: UserModel = Depends(get_active_user),
):
    """Update personal profile."""
    data = jsonable_encoder(user)
    params = UserUpdateDto(**data)
    if password is not None:
        params.password = password
    if name is not None:
        params.name = name
    if email is not None:
        params.email = email
    item = service.update(db, item=user, params=params)
    return item


@router.get("/users/profile", response_model=UserDto)
def get_me(db: Session = Depends(get_db), user: UserModel = Depends(get_active_user)):
    """Get current user."""
    return user


@router.get("/users/{id}", response_model=UserDto)
def get_by_id(
    id: int,
    db: Session = Depends(get_db),
    current_user: UserModel = Depends(get_active_user),
):
    """Get user by id."""
    user = service.get_by_id(db, id)
    if user == current_user:
        return user
    if not current_user.is_admin:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="The user doesn't have enough privileges")
    return user


@router.put("/users/{id}", response_model=UserDto)
def update(
    id: int,
    params: UserUpdateDto,
    db: Session = Depends(get_db),
    user: UserModel = Depends(get_admin),
):
    """Update user."""
    item = service.get_by_id(db, id)
    if not item:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="The user with this username does not exist in the system",
        )
    item = service.update(db, item=item, params=params)
    return item


@router.get("/users/check/{email}")
def check_user_exists(email: str, db: Session = Depends(get_db)):
    """Check if user with the email exists."""
    user = service.get_by_email(db, email=email)
    if user:
        return {"exists": True}
    else:
        return {"exists": False}
