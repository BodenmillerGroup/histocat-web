import copy
import logging
from http import HTTPStatus
from typing import Sequence

from fastapi import HTTPException

from histocat.cellxgene.common.app_config import AppConfig
from histocat.cellxgene.common.constants import Axis
from histocat.cellxgene.common.data_locator import DataLocator
from histocat.cellxgene.data_anndata.anndata_adaptor import AnndataAdaptor

logger = logging.getLogger(__name__)


def get_schema(app_config: AppConfig):
    """helper function to gather the schema from the data source and annotations"""
    data_locator = DataLocator("/app/histocat/data/pbmc3k.h5ad")
    data_adaptor = AnndataAdaptor(data_locator, app_config)
    annotations = None

    schema = data_adaptor.get_schema()
    schema = copy.deepcopy(schema)

    # add label obs annotations as needed
    if annotations is not None:
        label_schema = annotations.get_schema(data_adaptor)
        schema["annotations"]["obs"]["columns"].extend(label_schema)

    return schema


def get_annotations_obs(app_config: AppConfig, fields: Sequence[str]):
    data_locator = DataLocator("/app/histocat/data/pbmc3k.h5ad")
    data_adaptor = AnndataAdaptor(data_locator, app_config)
    annotations = None

    num_columns_requested = len(data_adaptor.get_obs_keys()) if len(fields) == 0 else len(fields)
    if data_adaptor.config.exceeds_limit("column_request_max", num_columns_requested):
        raise HTTPException(status_code=HTTPStatus.BAD_REQUEST, detail=f"column_request_max limit exceed.")

    try:
        labels = None
        if annotations:
            labels = annotations.read_labels(data_adaptor)
        fbs = data_adaptor.annotation_to_fbs_matrix(Axis.OBS, fields, labels)
        return fbs
    except KeyError as e:
        raise HTTPException(status_code=HTTPStatus.BAD_REQUEST, detail=e)


def get_annotations_var(app_config: AppConfig, fields: Sequence[str]):
    data_locator = DataLocator("/app/histocat/data/pbmc3k.h5ad")
    data_adaptor = AnndataAdaptor(data_locator, app_config)
    annotations = None

    num_columns_requested = len(data_adaptor.get_obs_keys()) if len(fields) == 0 else len(fields)
    if data_adaptor.config.exceeds_limit("column_request_max", num_columns_requested):
        raise HTTPException(status_code=HTTPStatus.BAD_REQUEST, detail=f"column_request_max limit exceed.")

    try:
        labels = None
        if annotations:
            labels = annotations.read_labels(data_adaptor)
        fbs = data_adaptor.annotation_to_fbs_matrix(Axis.VAR, fields, labels)
        return fbs
    except KeyError as e:
        raise HTTPException(status_code=HTTPStatus.BAD_REQUEST, detail=e)
