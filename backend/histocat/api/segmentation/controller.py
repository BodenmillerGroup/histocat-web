import logging

import dramatiq
from fastapi import APIRouter, Depends
from fastapi.responses import ORJSONResponse
from sqlalchemy.orm import Session

from histocat.api.db import get_db
from histocat.api.security import get_active_member
from histocat.core.member.models import MemberModel
from histocat.core.segmentation.dto import SegmentationSubmissionDto

logger = logging.getLogger(__name__)

router = APIRouter()


@router.post("/groups/{group_id}/projects/{project_id}/segmentation/process")
def process_segmentation(
    group_id: int,
    project_id: int,
    params: SegmentationSubmissionDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """
    Start segmentation processing
    """
    broker = dramatiq.get_broker()
    message = dramatiq.Message(
        actor_name="process_segmentation",
        queue_name="process",
        args=(),
        kwargs={"project_id": project_id, "json": params.json(),},
        options={},
    )
    broker.enqueue(message)
    return ORJSONResponse({"status": "submitted"})
