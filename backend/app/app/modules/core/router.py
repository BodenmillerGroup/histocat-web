from fastapi import APIRouter, Depends
from pydantic.types import EmailStr

import app.worker as worker
from app.api.utils.security import get_current_active_superuser
from app.core.utils import send_test_email
from app.modules.user.models import UserDBModel
from .models import MsgModel

router = APIRouter()


@router.post("/test-worker/", response_model=MsgModel, status_code=201)
def test_worker(
    msg: MsgModel, current_user: UserDBModel = Depends(get_current_active_superuser)
):
    """
    Test worker
    """
    worker.test_worker.send(msg.msg)
    return {"msg": "Word received"}


@router.post("/test-email/", response_model=MsgModel, status_code=201)
def test_email(
    email_to: EmailStr,
    current_user: UserDBModel = Depends(get_current_active_superuser),
):
    """
    Test emails
    """
    send_test_email(email_to=email_to)
    return {"msg": "Test email submitted"}
