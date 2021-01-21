from typing import Sequence

import dramatiq
from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import ORJSONResponse
from sqlalchemy.orm import Session

from histocat.api.db import get_db
from histocat.api.security import get_active_member, get_active_user
from histocat.core.member.models import MemberModel
from histocat.core.pipeline import service
from histocat.core.pipeline.dto import (
    PipelineCreateDto,
    PipelineDto,
    PipelineProcessDto,
    PipelineUpdateDto,
)
from histocat.core.user.models import UserModel

router = APIRouter()


@router.get("/pipelines/{pipeline_id}", response_model=PipelineDto)
def get_by_id(
    pipeline_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Get pipeline by id
    """
    item = service.get_by_id(db, id=pipeline_id)
    return item


@router.get("/projects/{project_id}/pipelines", response_model=Sequence[PipelineDto])
def get_project_pipelines(
    project_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Get all project pipelines
    """
    item = service.get_project_pipelines(db, project_id=project_id)
    return item


@router.post("/groups/{group_id}/pipelines", response_model=PipelineDto)
def create(
    group_id: int,
    params: PipelineCreateDto,
    db: Session = Depends(get_db),
    member: MemberModel = Depends(get_active_member),
):
    """
    Create new pipeline
    """
    items = service.create(db, params=params)
    return items


@router.patch("/groups/{group_id}/pipelines/{pipeline_id}", response_model=PipelineDto)
def update(
    group_id: int,
    pipeline_id: int,
    params: PipelineUpdateDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """
    Update pipeline
    """
    item = service.get_by_id(db, id=pipeline_id)
    if not item:
        raise HTTPException(
            status_code=404, detail="The pipeline with this id does not exist.",
        )
    item = service.update(db, item=item, params=params)
    return item


@router.delete("/pipelines/{pipeline_id}", response_model=int)
def delete_by_id(
    pipeline_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Delete pipeline by id
    """
    service.delete_by_id(db, id=pipeline_id)
    return pipeline_id


@router.post("/groups/{group_id}/pipelines/process")
def process(
    group_id: int,
    params: PipelineProcessDto,
    db: Session = Depends(get_db),
    member: MemberModel = Depends(get_active_member),
):
    """
    Process pipeline
    """
    broker = dramatiq.get_broker()
    message = dramatiq.Message(
        actor_name="process_pipeline",
        queue_name="process",
        args=(),
        kwargs={"dataset_id": params.dataset_id, "acquisition_ids": params.acquisition_ids, "steps": params.steps},
        options={},
    )
    broker.enqueue(message)
    return ORJSONResponse({"status": "submitted"})
