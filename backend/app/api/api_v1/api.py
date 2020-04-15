from fastapi import APIRouter

from app.modules.acquisition import controller as acquisition_controller
from app.modules.analysis import controller as analysis_controller
from app.modules.auth import controller as auth_controller
from app.modules.core import controller as core_controller
from app.modules.dataset import controller as dataset_controller
from app.modules.experiment import controller as experiment_controller
from app.modules.panorama import controller as panorama_controller
from app.modules.share import controller as share_controller
from app.modules.slide import controller as slide_controller
from app.modules.user import controller as user_controller

api_router = APIRouter()
api_router.include_router(auth_controller.router, prefix="/auth", tags=["auth"])
api_router.include_router(core_controller.router, prefix="/utils", tags=["utils"])
api_router.include_router(experiment_controller.router, prefix="/experiments", tags=["experiments"])
api_router.include_router(slide_controller.router, prefix="/slides", tags=["slides"])
api_router.include_router(panorama_controller.router, prefix="/panoramas", tags=["panoramas"])
api_router.include_router(acquisition_controller.router, prefix="/acquisitions", tags=["acquisitions"])
api_router.include_router(user_controller.router, prefix="/users", tags=["users"])
api_router.include_router(dataset_controller.router, prefix="/datasets", tags=["datasets"])
api_router.include_router(share_controller.router, prefix="/share", tags=["share"])
api_router.include_router(analysis_controller.router, prefix="/analysis", tags=["analysis"])
