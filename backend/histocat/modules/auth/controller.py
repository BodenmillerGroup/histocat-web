from datetime import timedelta

from fastapi import APIRouter, Body, Depends, HTTPException
from fastapi.security import OAuth2PasswordRequestForm
from pydantic import EmailStr
from sqlalchemy.orm import Session
from starlette import status
from starlette.responses import RedirectResponse

from histocat.api.db import get_db
from histocat.api.security import get_user
from histocat.config import config
from histocat.core.jwt import create_access_token
from histocat.core.security import get_password_hash
from histocat.core.utils import (
    generate_confirm_signup_token,
    generate_password_reset_token,
    send_confirm_signup_email,
    send_reset_password_email,
    verify_confirm_signup_token,
    verify_password_reset_token,
)
from histocat.modules.core.dto import MsgDto
from histocat.modules.user import service
from histocat.modules.user.dto import UserCreateDto, UserDto
from histocat.modules.user.models import UserModel as DBUser

from .dto import TokenDto

router = APIRouter()


@router.post("/auth/login", response_model=TokenDto)
def login_access_token(db: Session = Depends(get_db), form_data: OAuth2PasswordRequestForm = Depends()):
    """
    OAuth2 compatible token login, get an access token for future requests
    """
    user = service.authenticate(db, email=form_data.username, password=form_data.password)
    if not user:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Incorrect email or password")
    elif not user.is_active:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Inactive user")
    access_token_expires = timedelta(minutes=config.ACCESS_TOKEN_EXPIRE_MINUTES)
    return {
        "access_token": create_access_token(data={"user_id": user.id}, expires_delta=access_token_expires),
        "token_type": "bearer",
    }


@router.post("/auth/test-token", response_model=UserDto)
def test_token(user: DBUser = Depends(get_user)):
    """
    Test access token
    """
    return user


@router.post("/auth/password-recovery/{email}", response_model=MsgDto)
def recover_password(email: str, db: Session = Depends(get_db)):
    """
    Password Recovery
    """
    user = service.get_by_email(db, email=email)

    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="The user with this username does not exist in the system.",
        )
    password_reset_token = generate_password_reset_token(email=email)
    send_reset_password_email(email_to=user.email, email=email, token=password_reset_token)
    return {"msg": "Password recovery email sent"}


@router.post("/auth/reset-password/", response_model=MsgDto)
def reset_password(token: str = Body(...), new_password: str = Body(...), db: Session = Depends(get_db)):
    """
    Reset password
    """
    email = verify_password_reset_token(token)
    if not email:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid token")
    user = service.get_by_email(db, email=email)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="The user with this username does not exist in the system.",
        )
    elif not user.is_active:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Inactive user")
    hashed_password = get_password_hash(new_password)
    user.password = hashed_password
    db.add(user)
    db.commit()
    return {"msg": "Password updated successfully"}


@router.post("/auth/signup", response_model=UserDto)
def signup(
    password: str = Body(...), email: EmailStr = Body(...), name: str = Body(None), db: Session = Depends(get_db),
):
    """Create new user without the need to be logged in."""
    if not config.USERS_OPEN_REGISTRATION:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN, detail="Open user resgistration is forbidden on this server",
        )
    user = service.get_by_email(db, email=email)
    if user:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail="The user with this username already exists",
        )

    confirm_signup_token = generate_confirm_signup_token(email=email)
    send_confirm_signup_email(email_to=email, username=name, token=confirm_signup_token)

    params = UserCreateDto(password=password, email=email, name=name, is_active=False)
    service.create(db, params=params)


@router.get("/auth/confirm-signup/{token}")
def confirm_signup(token: str, db: Session = Depends(get_db)):
    email = verify_confirm_signup_token(token)
    if not email:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid token")
    user = service.get_by_email(db, email=email)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="User doesn't exist",
        )
    service.confirm_signup(db, item=user)
    return RedirectResponse("/")
