import dramatiq
from fastapi import APIRouter, Depends
from pydantic import EmailStr

from histocat.api.security import get_admin
from histocat.core.utils import send_test_email

from .dto import MsgDto

router = APIRouter()


@router.post("/test-worker/", response_model=MsgDto, status_code=201)
def test_worker(msg: MsgDto, user=Depends(get_admin)):
    """
    Test worker
    """
    broker = dramatiq.get_broker()
    message = dramatiq.Message(
        actor_name="test_worker", queue_name="default", args=(), kwargs={"word": msg.msg}, options={},
    )
    broker.enqueue(message)
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
