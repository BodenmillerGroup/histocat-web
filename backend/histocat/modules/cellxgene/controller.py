from fastapi import APIRouter

from . import service

router = APIRouter()


# Initialization routes

@router.get("/schema", status_code=200)
def get_schema():
    schema = service.get_schema()
    return { "schema": schema }


@router.get("/config", status_code=200)
def get_config():
    return {"result": "config"}


# Data routes

@router.get("/annotations/obs", status_code=200)
def get_annotations_obs():
    return {"result": "annotations_obs"}


@router.get("/annotations/var", status_code=200)
def get_annotations_var():
    return {"result": "annotations_var"}


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
