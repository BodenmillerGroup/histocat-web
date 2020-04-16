from fastapi import APIRouter, Depends
from fastapi.responses import ORJSONResponse
from pydantic import EmailStr

import app.worker as worker
from app.api.utils.security import get_current_active_superuser
from app.core.utils import send_test_email

from .dto import MsgDto

router = APIRouter()


@router.post("/test-worker/", response_model=MsgDto, status_code=201)
def test_worker(msg: MsgDto, current_user=Depends(get_current_active_superuser)):
    """
    Test worker
    """
    worker.test_worker.send(msg.msg)
    return ORJSONResponse({"msg": "Word received"})


@router.post("/test-email/", response_model=MsgDto, status_code=201)
def test_email(
    email_to: EmailStr, current_user=Depends(get_current_active_superuser),
):
    """
    Test emails
    """
    send_test_email(email_to=email_to)
    return ORJSONResponse({"msg": "Test email submitted"})
