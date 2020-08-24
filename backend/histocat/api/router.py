from fastapi import APIRouter

from histocat.modules.acquisition import controller as acquisition_controller
from histocat.modules.analysis import controller as analysis_controller
from histocat.modules.auth import controller as auth_controller
from histocat.modules.core import controller as core_controller
from histocat.modules.dataset import controller as dataset_controller
from histocat.modules.experiment import controller as experiment_controller
from histocat.modules.gate import controller as gate_controller
from histocat.modules.group import controller as group_controller
from histocat.modules.member import controller as member_controller
from histocat.modules.panorama import controller as panorama_controller
from histocat.modules.preset import controller as preset_controller
from histocat.modules.slide import controller as slide_controller
from histocat.modules.user import controller as user_controller

api_router = APIRouter()

api_router.include_router(auth_controller.router, tags=["auth"])
api_router.include_router(core_controller.router, tags=["utils"])
api_router.include_router(group_controller.router, tags=["groups"])
api_router.include_router(experiment_controller.router, tags=["experiments"])
api_router.include_router(slide_controller.router, tags=["slides"])
api_router.include_router(panorama_controller.router, tags=["panoramas"])
api_router.include_router(acquisition_controller.router, tags=["acquisitions"])
api_router.include_router(user_controller.router, tags=["users"])
api_router.include_router(member_controller.router, tags=["members"])
api_router.include_router(dataset_controller.router, tags=["datasets"])
api_router.include_router(preset_controller.router, tags=["presets"])
api_router.include_router(gate_controller.router, tags=["gates"])
api_router.include_router(analysis_controller.router, tags=["analysis"])
