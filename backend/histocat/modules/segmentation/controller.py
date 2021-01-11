import logging

from fastapi import APIRouter, Depends
from fastapi.responses import ORJSONResponse
from sqlalchemy.orm import Session

from histocat import worker
from histocat.api.db import get_db
from histocat.api.security import get_active_member
from histocat.modules.member.models import MemberModel

from .dto import SegmentationSubmissionDto

logger = logging.getLogger(__name__)

router = APIRouter()


@router.post("/groups/{group_id}/projects/{project_id}/segmentation/process")
def process_segmentation(
    group_id: int,
    project_id: int,
    params: SegmentationSubmissionDto, member: MemberModel = Depends(get_active_member), db: Session = Depends(get_db),
):
    """
    Start segmentation processing
    """

    worker.process_segmentation.send(group_id, project_id, params.dict())
    return ORJSONResponse({"status": "submitted"})
