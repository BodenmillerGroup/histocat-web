import jwt
from fastapi import Depends, HTTPException, Security
from fastapi.security import OAuth2PasswordBearer
from jwt import PyJWTError
from sqlalchemy.orm import Session
from starlette.status import HTTP_403_FORBIDDEN

from histocat.api.utils.db import get_db
from histocat.core import config
from histocat.core.jwt import ALGORITHM
from histocat.modules.auth.dto import TokenPayloadDto
from histocat.modules.user import service
from histocat.modules.user.models import UserModel

reusable_oauth2 = OAuth2PasswordBearer(tokenUrl="/api/v1/auth/login")


def get_current_user(db: Session = Depends(get_db), token: str = Security(reusable_oauth2)):
    try:
        payload = jwt.decode(token, config.SECRET_KEY, algorithms=[ALGORITHM])
        token_data = TokenPayloadDto(**payload)
    except PyJWTError:
        raise HTTPException(status_code=HTTP_403_FORBIDDEN, detail="Could not validate credentials")
    user = service.get_by_id(db, token_data.user_id)
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    return user


def get_current_active_user(user: UserModel = Security(get_current_user)):
    if not user.is_active:
        raise HTTPException(status_code=400, detail="Inactive user")
    return user


def get_current_active_superuser(user: UserModel = Security(get_current_user)):
    if not user.is_admin:
        raise HTTPException(status_code=400, detail="The user doesn't have enough privileges")
    return user
