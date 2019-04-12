from typing import List

from fastapi import APIRouter, Body, Depends, HTTPException
from fastapi.encoders import jsonable_encoder
from pydantic.types import EmailStr
from sqlalchemy.orm import Session

from app import crud
from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.core import config
from app.db_models.user import User
from app.models.user import UserModel, UserInCreateModel, UserInDBModel, UserInUpdateModel
from app.utils import send_new_account_email

router = APIRouter()


@router.get("/users/", tags=["users"], response_model=List[UserModel])
def read_users(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Retrieve users
    """

    users = crud.user.get_multi(db, skip=skip, limit=limit)
    return users


@router.post("/users/", tags=["users"], response_model=UserModel)
def create_user(
    *,
    db: Session = Depends(get_db),
    params: UserInCreateModel,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Create new user
    """
    user = crud.user.get_by_email(db, email=params.email)
    if user:
        raise HTTPException(
            status_code=400,
            detail="The user with this username already exists in the system.",
        )
    user = crud.user.create(db, params=params)
    if config.EMAILS_ENABLED and params.email:
        send_new_account_email(
            email_to=params.email, username=params.email, password=params.password
        )
    return user


@router.put("/users/me", tags=["users"], response_model=UserModel)
def update_user_me(
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
    current_user_data = jsonable_encoder(current_user)
    user_in = UserInUpdateModel(**current_user_data)
    if password is not None:
        user_in.password = password
    if full_name is not None:
        user_in.full_name = full_name
    if email is not None:
        user_in.email = email
    user = crud.user.update(db, user=current_user, params=user_in)
    return user


@router.get("/users/me", tags=["users"], response_model=UserModel)
def read_user_me(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user),
):
    """
    Get current user
    """
    return current_user


@router.post("/users/open", tags=["users"], response_model=UserModel)
def create_user_open(
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
            status_code=403,
            detail="Open user resgistration is forbidden on this server",
        )
    user = crud.user.get_by_email(db, email=email)
    if user:
        raise HTTPException(
            status_code=400,
            detail="The user with this username already exists in the system",
        )
    user_in = UserInCreateModel(password=password, email=email, full_name=full_name)
    user = crud.user.create(db, params=user_in)
    return user


@router.get("/users/{user_id}", tags=["users"], response_model=UserModel)
def read_user_by_id(
    user_id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get a specific user by id
    """
    user = crud.user.get(db, id=user_id)
    if user == current_user:
        return user
    if not crud.user.is_superuser(current_user):
        raise HTTPException(
            status_code=400, detail="The user doesn't have enough privileges"
        )
    return user


@router.put("/users/{user_id}", tags=["users"], response_model=UserModel)
def update_user(
    *,
    db: Session = Depends(get_db),
    user_id: int,
    user_in: UserInUpdateModel,
    current_user: UserInDBModel = Depends(get_current_active_superuser),
):
    """
    Update a user
    """
    user = crud.user.get(db, id=user_id)

    if not user:
        raise HTTPException(
            status_code=404,
            detail="The user with this username does not exist in the system",
        )
    user = crud.user.update(db, user=user, params=user_in)
    return user
