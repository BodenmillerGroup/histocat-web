from fastapi import APIRouter

from histocat.api.acquisition import controller as acquisition_controller
from histocat.api.analysis import controller as analysis_controller
from histocat.api.auth import controller as auth_controller
from histocat.api.check import controller as check_controller
from histocat.api.dataset import controller as dataset_controller
from histocat.api.gate import controller as gate_controller
from histocat.api.group import controller as group_controller
from histocat.api.member import controller as member_controller
from histocat.api.model import controller as model_controller
from histocat.api.panorama import controller as panorama_controller
from histocat.api.pipeline import controller as pipeline_controller
from histocat.api.preset import controller as preset_controller
from histocat.api.project import controller as project_controller
from histocat.api.result import controller as result_controller
from histocat.api.segmentation import controller as segmentation_controller
from histocat.api.slide import controller as slide_controller
from histocat.api.user import controller as user_controller

api_router = APIRouter()

api_router.include_router(auth_controller.router, tags=["auth"])
api_router.include_router(check_controller.router, tags=["check"])
api_router.include_router(group_controller.router, tags=["groups"])
api_router.include_router(project_controller.router, tags=["projects"])
api_router.include_router(slide_controller.router, tags=["slides"])
api_router.include_router(panorama_controller.router, tags=["panoramas"])
api_router.include_router(acquisition_controller.router, tags=["acquisitions"])
api_router.include_router(user_controller.router, tags=["users"])
api_router.include_router(member_controller.router, tags=["members"])
api_router.include_router(dataset_controller.router, tags=["datasets"])
api_router.include_router(pipeline_controller.router, tags=["pipelines"])
api_router.include_router(result_controller.router, tags=["results"])
api_router.include_router(preset_controller.router, tags=["presets"])
api_router.include_router(gate_controller.router, tags=["gates"])
api_router.include_router(analysis_controller.router, tags=["analysis"])
api_router.include_router(model_controller.router, tags=["models"])
api_router.include_router(segmentation_controller.router, tags=["segmentation"])
