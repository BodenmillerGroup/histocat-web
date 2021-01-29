import logging
import os
import uuid
from typing import Sequence

import dramatiq
from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile
from sqlalchemy.orm import Session
from starlette import status

from histocat.api.db import get_db
from histocat.api.security import get_active_member, get_group_admin
from histocat.config import config
from histocat.core.member.models import MemberModel
from histocat.core.model import service
from histocat.core.model.dto import ModelCreateDto, ModelDto, ModelUpdateDto

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/groups/{group_id}/models", response_model=Sequence[ModelDto])
def get_group_models(group_id: int, db: Session = Depends(get_db), member: MemberModel = Depends(get_active_member)):
    """
    Get group models
    """
    items = service.get_group_models(db, group_id=group_id)
    return items


@router.get("/groups/{group_id}/models/{model_id}", response_model=ModelDto)
def get_by_id(
    group_id: int, model_id: int, member: MemberModel = Depends(get_active_member), db: Session = Depends(get_db),
):
    """
    Get model by id
    """
    item = service.get(db, id=model_id)
    return item


@router.delete("/groups/{group_id}/models/{model_id}", response_model=int)
def delete_by_id(
    group_id: int, model_id: int, member: MemberModel = Depends(get_group_admin), db: Session = Depends(get_db),
):
    """
    Delete model by id
    """
    service.remove(db, id=model_id)
    return model_id


@router.patch("/groups/{group_id}/models/{model_id}", response_model=ModelDto)
def update(
    group_id: int,
    model_id: int,
    params: ModelUpdateDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """
    Update model
    """
    item = service.get(db, id=model_id)
    if not item:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Model id:{model_id} not found")
    item = service.update(db, item=item, params=params)
    return item


@router.post("/groups/{group_id}/models", response_model=ModelDto)
def create(
    group_id: int,
    name: str = Form(""),
    description: str = Form(None),
    file: UploadFile = File(None),
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    item = service.get_by_group_id_and_name(db, group_id=group_id, name=name)
    if item:
        raise HTTPException(
            status_code=400, detail="The model with this name already exists.",
        )

    params = ModelCreateDto(name=name, description=description)
    model = service.create(db, group_id=group_id, params=params)

    path = os.path.join(config.INBOX_DIRECTORY, str(uuid.uuid4()))
    if not os.path.exists(path):
        os.makedirs(path)
    uri = os.path.join(path, file.filename)
    with open(uri, "wb") as f:
        f.write(file.file.read())

    broker = dramatiq.get_broker()
    message = dramatiq.Message(
        actor_name="import_model",
        queue_name="import",
        args=(),
        kwargs={"uri": uri, "group_id": group_id, "model_id": model.id},
        options={},
    )
    broker.enqueue(message)
    return model
