from fastapi import APIRouter, Depends
from pydantic import EmailStr

import histocat.worker as worker
from histocat.api.security import get_admin
from histocat.core.utils import send_test_email

from .dto import MsgDto

router = APIRouter()


@router.post("/test-worker/", response_model=MsgDto, status_code=201)
def test_worker(msg: MsgDto, user=Depends(get_admin)):
    """
    Test worker
    """
    worker.test_worker.send(msg.msg)
    return {"msg": "Word received"}


@router.post("/test-email/", response_model=MsgDto, status_code=201)
def test_email(
    email_to: EmailStr, user=Depends(get_admin),
):
    """
    Test emails
    """
    send_test_email(email_to=email_to)
    return {"msg": "Test email submitted"}
