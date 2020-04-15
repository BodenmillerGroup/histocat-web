from typing import Sequence

from fastapi import APIRouter, Body, Depends, HTTPException
from fastapi.encoders import jsonable_encoder
from pydantic import EmailStr
from sqlalchemy.orm import Session

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.core import config
from app.core.utils import send_new_account_email
from app.modules.user.models import UserModel

from . import service
from .dto import UserCreateDto, UserDto, UserUpdateDto

router = APIRouter()


@router.get("/", response_model=Sequence[UserDto])
def get_all(db: Session = Depends(get_db), current_user: UserModel = Depends(get_current_active_user)):
    """Get all users."""
    items = service.get_all(db)
    return items


@router.post("/", response_model=UserDto)
def create(
    params: UserCreateDto,
    db: Session = Depends(get_db),
    current_user: UserModel = Depends(get_current_active_superuser),
):
    """Create new user."""
    item = service.get_by_email(db, email=params.email)
    if item:
        raise HTTPException(
            status_code=400, detail="The user with this username already exists in the system.",
        )
    item = service.create(db, params=params)
    if config.EMAILS_ENABLED and params.email:
        send_new_account_email(email_to=params.email, username=params.email, password=params.password)
    return item


@router.patch("/profile", response_model=UserDto)
def update_me(
    password: str = Body(None),
    name: str = Body(None),
    email: EmailStr = Body(None),
    db: Session = Depends(get_db),
    current_user: UserModel = Depends(get_current_active_user),
):
    """Update personal profile."""
    data = jsonable_encoder(current_user)
    params = UserUpdateDto(**data)
    if password is not None:
        params.password = password
    if name is not None:
        params.name = name
    if email is not None:
        params.email = email
    item = service.update(db, item=current_user, params=params)
    return item


@router.get("/profile", response_model=UserDto)
def get_me(db: Session = Depends(get_db), current_user: UserModel = Depends(get_current_active_user)):
    """Get current user."""
    return current_user


@router.post("/signup", response_model=UserDto)
def create_open(
    password: str = Body(...), email: EmailStr = Body(...), name: str = Body(None), db: Session = Depends(get_db),
):
    """Create new user without the need to be logged in."""
    if not config.USERS_OPEN_REGISTRATION:
        raise HTTPException(
            status_code=403, detail="Open user resgistration is forbidden on this server",
        )
    user = service.get_by_email(db, email=email)
    if user:
        raise HTTPException(
            status_code=400, detail="The user with this username already exists in the system",
        )
    user_in = UserCreateDto(password=password, email=email, name=name)
    item = service.create(db, params=user_in)
    return item


@router.get("/{id}", response_model=UserDto)
def get_by_id(
    id: int, db: Session = Depends(get_db), current_user: UserModel = Depends(get_current_active_user),
):
    """Get user by id."""
    user = service.get_by_id(db, id)
    if user == current_user:
        return user
    if not current_user.is_admin:
        raise HTTPException(status_code=400, detail="The user doesn't have enough privileges")
    return user


@router.put("/{id}", response_model=UserDto)
def update(
    id: int,
    params: UserUpdateDto,
    db: Session = Depends(get_db),
    current_user: UserModel = Depends(get_current_active_superuser),
):
    """Update user."""
    item = service.get_by_id(db, id)
    if not item:
        raise HTTPException(
            status_code=404, detail="The user with this username does not exist in the system",
        )
    item = service.update(db, item=item, params=params)
    return item


@router.get("/check/{email}")
def check_user_exists(email: str, db: Session = Depends(get_db)):
    """Check if user with the email exists."""
    user = service.get_by_email(db, email=email)
    if user:
        return {"exists": True}
    else:
        return {"exists": False}
