import jwt
from fastapi import Depends, HTTPException, Security
from fastapi.security import OAuth2PasswordBearer
from jwt import PyJWTError
from sqlalchemy.orm import Session
from starlette import status

from histocat.api.db import get_db
from histocat.config import config
from histocat.core.auth.dto import TokenPayloadDto
from histocat.core.jwt import ALGORITHM
from histocat.core.member import service as member_service
from histocat.core.member.models import MemberModel
from histocat.core.user import service as user_service
from histocat.core.user.models import UserModel

reusable_oauth2 = OAuth2PasswordBearer(tokenUrl="/api/v1/auth/login")


def get_token_data(token: str = Security(reusable_oauth2)):
    try:
        payload = jwt.decode(token, config.JWT_SECRET, algorithms=[ALGORITHM])
        token_data = TokenPayloadDto(**payload)
    except PyJWTError:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Could not validate credentials")
    return token_data


def get_user(db: Session = Depends(get_db), token_data: TokenPayloadDto = Depends(get_token_data)):
    user = user_service.get_by_id(db, token_data.user_id)
    if not user:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="User not found")
    return user


def get_active_user(user: UserModel = Security(get_user)):
    if not user.is_active:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="User account is deactivated")
    return user


def get_admin(user: UserModel = Security(get_user)):
    if not user.is_active or not user.is_admin:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="User doesn't have enough privileges")
    return user


def get_group_member(
    group_id: int, db: Session = Depends(get_db), token_data: TokenPayloadDto = Depends(get_token_data)
):
    member = member_service.get_by_group_id_and_user_id(db, group_id=group_id, user_id=token_data.user_id)
    if not member:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Group member not found")
    if not member.user.is_active:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="User account is deactivated")
    return member


def get_active_member(member: MemberModel = Security(get_group_member)):
    if not member.is_active:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Group member account is deactivated")
    return member


def get_group_admin(member: MemberModel = Security(get_group_member)):
    if not member.is_active or member.role < 100:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED, detail="Group member doesn't have enough privileges"
        )
    return member
