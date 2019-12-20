from typing import List

from fastapi import APIRouter, Body, Depends, HTTPException
from fastapi.encoders import jsonable_encoder
from pydantic import EmailStr
from sqlalchemy.orm import Session

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.core import config
from app.core.utils import send_new_account_email
from app.modules.user.db import User

from . import crud
from .models import UserCreateModel, UserDBModel, UserModel, UserUpdateModel

router = APIRouter()


@router.get("/", response_model=List[UserModel])
def read_all(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_user),
):
    """
    Retrieve users
    """
    items = crud.get_multi(db, skip=skip, limit=limit)
    return items


@router.post("/", response_model=UserModel)
def create(
    *,
    db: Session = Depends(get_db),
    params: UserCreateModel,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Create new user
    """
    item = crud.get_by_email(db, email=params.email)
    if item:
        raise HTTPException(
            status_code=400, detail="The user with this username already exists in the system.",
        )
    item = crud.create(db, params=params)
    if config.EMAILS_ENABLED and params.email:
        send_new_account_email(email_to=params.email, username=params.email, password=params.password)
    return item


@router.patch("/profile", response_model=UserModel)
def update_me(
    *,
    db: Session = Depends(get_db),
    password: str = Body(None),
    full_name: str = Body(None),
    email: EmailStr = Body(None),
    current_user: User = Depends(get_current_active_user),
):
    """
    Update own user
    """
    data = jsonable_encoder(current_user)
    params = UserUpdateModel(**data)
    if password is not None:
        params.password = password
    if full_name is not None:
        params.full_name = full_name
    if email is not None:
        params.email = email
    item = crud.update(db, item=current_user, params=params)
    return item


@router.get("/profile", response_model=UserModel)
def read_me(db: Session = Depends(get_db), current_user: User = Depends(get_current_active_user)):
    """
    Get current user
    """
    return current_user


@router.post("/signup", response_model=UserModel)
def create_open(
    *,
    db: Session = Depends(get_db),
    password: str = Body(...),
    email: EmailStr = Body(...),
    full_name: str = Body(None),
):
    """
    Create new user without the need to be logged in
    """
    if not config.USERS_OPEN_REGISTRATION:
        raise HTTPException(
            status_code=403, detail="Open user resgistration is forbidden on this server",
        )
    user = crud.get_by_email(db, email=email)
    if user:
        raise HTTPException(
            status_code=400, detail="The user with this username already exists in the system",
        )
    user_in = UserCreateModel(password=password, email=email, full_name=full_name)
    item = crud.create(db, params=user_in)
    return item


@router.get("/{id}", response_model=UserModel)
def read_by_id(
    id: int, current_user: User = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Get a specific user by id
    """
    user = crud.get(db, id=id)
    if user == current_user:
        return user
    if not crud.is_superuser(current_user):
        raise HTTPException(status_code=400, detail="The user doesn't have enough privileges")
    return user


@router.put("/{id}", response_model=UserModel)
def update(
    *,
    db: Session = Depends(get_db),
    id: int,
    params: UserUpdateModel,
    current_user: UserDBModel = Depends(get_current_active_superuser),
):
    """
    Update a user
    """
    item = crud.get(db, id=id)

    if not item:
        raise HTTPException(
            status_code=404, detail="The user with this username does not exist in the system",
        )
    item = crud.update(db, item=item, params=params)
    return item


@router.get("/check/{email}")
def check_user_exists(email: str, db: Session = Depends(get_db)):
    """
    Check if user with the email exists
    """
    user = crud.get_by_email(db, email=email)
    if user:
        return {"exists": True}
    else:
        return {"exists": False}
