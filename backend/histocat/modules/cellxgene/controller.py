from typing import Sequence

from fastapi import APIRouter, Query
from starlette.requests import Request
from starlette.responses import StreamingResponse

from . import service
from ...core.utils import stream_bytes

router = APIRouter()


# Initialization routes


@router.get("/schema", status_code=200)
async def get_schema(request: Request):
    schema = service.get_schema(request.app.app_config)
    return {"schema": schema}


@router.get("/config", status_code=200)
def get_config():
    return {"result": "config"}


# Data routes


@router.get("/annotations/obs", status_code=200)
def get_annotations_obs(request: Request, annotation_name: Sequence[str] = Query(None)):
    data = service.get_annotations_obs(request.app.app_config, annotation_name)
    return StreamingResponse(
        stream_bytes(data), media_type="application/octet-stream", headers={"content-length": str(len(data))}
    )


@router.get("/annotations/var", status_code=200)
def get_annotations_var(request: Request, annotation_name: Sequence[str] = Query(None)):
    data = service.get_annotations_var(request.app.app_config, annotation_name)
    return StreamingResponse(
        stream_bytes(data), media_type="application/octet-stream", headers={"content-length": str(len(data))}
    )


@router.get("/data/var", status_code=200)
def get_data_var():
    return {"result": "data_var"}


# Computation routes


@router.get("/diffexp/obs", status_code=200)
def get_diffexp_obs():
    return {"result": "diffexp_obs"}


@router.get("/layout/obs", status_code=200)
def get_layout_obs():
    return {"result": "layout_obs"}
