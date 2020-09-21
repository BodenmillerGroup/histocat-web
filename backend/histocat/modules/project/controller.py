import logging
import os
import uuid
from typing import Sequence, Set

from fastapi import APIRouter, Depends, File, HTTPException, UploadFile
from sqlalchemy.orm import Session

import histocat.worker as worker
from histocat.api.db import get_db
from histocat.api.security import get_active_member, get_group_admin
from histocat.config import config
from histocat.modules.member.models import MemberModel

from . import service
from .dto import ProjectCreateDto, ProjectDto, ProjectFullDto, ProjectUpdateDto

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/groups/{group_id}/tags", response_model=Set[str])
def get_tags(group_id: int, db: Session = Depends(get_db), member: MemberModel = Depends(get_active_member)):
    """
    Get projects tags
    """
    items = service.get_tags(db, group_id=group_id)
    return items


@router.get("/groups/{group_id}/projects", response_model=Sequence[ProjectDto])
def get_group_projects(group_id: int, db: Session = Depends(get_db), member: MemberModel = Depends(get_active_member)):
    """
    Get group projects
    """
    items = service.get_group_projects(db, group_id=group_id)
    return items


@router.post("/groups/{group_id}/projects", response_model=ProjectDto)
def create(
    group_id: int,
    params: ProjectCreateDto,
    db: Session = Depends(get_db),
    member: MemberModel = Depends(get_active_member),
):
    """
    Create new project
    """
    item = service.get_by_name(db, name=params.name)
    if item:
        raise HTTPException(
            status_code=400, detail="The project with this name already exists.",
        )
    item = service.create(db, group_id=group_id, params=params, member=member)
    return item


@router.get("/groups/{group_id}/projects/{project_id}", response_model=ProjectDto)
def get_by_id(
    group_id: int, project_id: int, member: MemberModel = Depends(get_active_member), db: Session = Depends(get_db),
):
    """
    Get project by id
    """
    item = service.get(db, id=project_id)
    return item


@router.delete("/groups/{group_id}/projects/{project_id}", response_model=ProjectDto)
def delete_by_id(
    group_id: int, project_id: int, member: MemberModel = Depends(get_group_admin), db: Session = Depends(get_db),
):
    """
    Delete project by id
    """
    item = service.remove(db, id=project_id)
    return item


@router.put("/groups/{group_id}/projects/{project_id}", response_model=ProjectDto)
def update(
    group_id: int,
    project_id: int,
    params: ProjectUpdateDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """
    Update project
    """
    item = service.get(db, id=project_id)
    if not item:
        raise HTTPException(
            status_code=404, detail="The project with this id does not exist.",
        )
    item = service.update(db, item=item, params=params)
    return item


@router.post("/groups/{group_id}/projects/{project_id}/upload")
def upload_data(
    group_id: int,
    project_id: int,
    file: UploadFile = File(None),
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    path = os.path.join(config.INBOX_DIRECTORY, str(uuid.uuid4()))
    if not os.path.exists(path):
        os.makedirs(path)
    uri = os.path.join(path, file.filename)
    with open(uri, "wb") as f:
        f.write(file.file.read())
    worker.import_data.send(uri, project_id)
    return {"uri": uri}


@router.get("/groups/{group_id}/projects/{project_id}/data", response_model=ProjectFullDto)
async def read_data(
    group_id: int, project_id: int, member: MemberModel = Depends(get_active_member), db: Session = Depends(get_db),
):
    """
    Get all project data
    """
    item = service.get_data(db, id=project_id)
    return item
