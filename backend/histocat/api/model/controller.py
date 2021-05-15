import logging
import os
import uuid
from typing import Sequence

import dramatiq
from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile
from sqlalchemy.orm import Session
from starlette import status

from histocat.api.db import get_db
from histocat.api.security import get_active_user, get_admin
from histocat.config import config
from histocat.core.model import service
from histocat.core.model.dto import ModelCreateDto, ModelDto, ModelUpdateDto
from histocat.core.user.models import UserModel

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/models", response_model=Sequence[ModelDto])
def get_all_models(db: Session = Depends(get_db), user: UserModel = Depends(get_active_user)):
    """
    Get all models
    """
    items = service.get_all_models(db)
    return items


@router.get("/models/{model_id}", response_model=ModelDto)
def get_by_id(
    model_id: int,
    user: UserModel = Depends(get_active_user),
    db: Session = Depends(get_db),
):
    """
    Get model by id
    """
    item = service.get(db, id=model_id)
    return item


@router.delete("/models/{model_id}", response_model=int)
def delete_by_id(
    model_id: int,
    user: UserModel = Depends(get_admin),
    db: Session = Depends(get_db),
):
    """
    Delete model by id
    """
    service.remove(db, id=model_id)
    return model_id


@router.patch("/models/{model_id}", response_model=ModelDto)
def update(
    model_id: int,
    params: ModelUpdateDto,
    user: UserModel = Depends(get_admin),
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


@router.post("/models", response_model=ModelDto)
def create(
    name: str = Form(""),
    application: str = Form(""),
    description: str = Form(None),
    file: UploadFile = File(None),
    user: UserModel = Depends(get_admin),
    db: Session = Depends(get_db),
):
    item = service.get_by_name(db, name=name)
    if item:
        raise HTTPException(
            status_code=400,
            detail="The model with this name already exists.",
        )

    params = ModelCreateDto(name=name, description=description, application=application)
    model = service.create(db, params=params)

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
        kwargs={"uri": uri, "model_id": model.id},
        options={},
    )
    broker.enqueue(message)
    return model
