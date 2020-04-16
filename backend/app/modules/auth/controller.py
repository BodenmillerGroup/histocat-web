from datetime import timedelta

from fastapi import APIRouter, Body, Depends, HTTPException
from fastapi.responses import ORJSONResponse
from fastapi.security import OAuth2PasswordRequestForm
from sqlalchemy.orm import Session

from app.api.utils.db import get_db
from app.api.utils.security import get_current_user
from app.core import config
from app.core.jwt import create_access_token
from app.core.security import get_password_hash
from app.core.utils import (
    generate_password_reset_token,
    send_reset_password_email,
    verify_password_reset_token,
)
from app.modules.core.dto import MsgDto
from app.modules.user import service
from app.modules.user.dto import UserDto
from app.modules.user.models import UserModel as DBUser

from .dto import TokenDto

router = APIRouter()


@router.post("/login", response_model=TokenDto)
def login_access_token(db: Session = Depends(get_db), form_data: OAuth2PasswordRequestForm = Depends()):
    """
    OAuth2 compatible token login, get an access token for future requests
    """
    user = service.authenticate(db, email=form_data.username, password=form_data.password)
    if not user:
        raise HTTPException(status_code=400, detail="Incorrect email or password")
    elif not user.is_active:
        raise HTTPException(status_code=400, detail="Inactive user")
    access_token_expires = timedelta(minutes=config.ACCESS_TOKEN_EXPIRE_MINUTES)
    return ORJSONResponse({
        "access_token": create_access_token(data={"user_id": user.id}, expires_delta=access_token_expires),
        "token_type": "bearer",
    })


@router.post("/test-token", response_model=UserDto)
def test_token(current_user: DBUser = Depends(get_current_user)):
    """
    Test access token
    """
    return ORJSONResponse(current_user)


@router.post("/password-recovery/{email}", response_model=MsgDto)
def recover_password(email: str, db: Session = Depends(get_db)):
    """
    Password Recovery
    """
    user = service.get_by_email(db, email=email)

    if not user:
        raise HTTPException(
            status_code=404, detail="The user with this username does not exist in the system.",
        )
    password_reset_token = generate_password_reset_token(email=email)
    send_reset_password_email(email_to=user.email, email=email, token=password_reset_token)
    return ORJSONResponse({"msg": "Password recovery email sent"})


@router.post("/reset-password/", response_model=MsgDto)
def reset_password(token: str = Body(...), new_password: str = Body(...), db: Session = Depends(get_db)):
    """
    Reset password
    """
    email = verify_password_reset_token(token)
    if not email:
        raise HTTPException(status_code=400, detail="Invalid token")
    user = service.get_by_email(db, email=email)
    if not user:
        raise HTTPException(
            status_code=404, detail="The user with this username does not exist in the system.",
        )
    elif not user.is_active:
        raise HTTPException(status_code=400, detail="Inactive user")
    hashed_password = get_password_hash(new_password)
    user.password = hashed_password
    db.add(user)
    db.commit()
    return ORJSONResponse({"msg": "Password updated successfully"})