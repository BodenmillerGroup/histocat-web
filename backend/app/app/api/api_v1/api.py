from fastapi import APIRouter

from app.modules.acquisition import router as acquisition_router
from app.modules.auth import router as auth_router
from app.modules.channel import router as channel_router
from app.modules.core import router as core_router
from app.modules.experiment import router as experiment_router
from app.modules.slide import router as slide_router
from app.modules.panorama import router as panorama_router
from app.modules.roi import router as roi_router
from app.modules.roi_point import router as roi_point_router
from app.modules.user import router as user_router
from app.modules.dataset import router as dataset_router
from app.modules.share import router as share_router
from app.modules.analysis import router as analysis_router

api_router = APIRouter()
api_router.include_router(auth_router.router, prefix="/auth", tags=["auth"])
api_router.include_router(core_router.router, prefix="/utils", tags=["utils"])
api_router.include_router(experiment_router.router, prefix="/experiments", tags=["experiments"])
api_router.include_router(slide_router.router, prefix="/slides", tags=["slides"])
api_router.include_router(panorama_router.router, prefix="/panoramas", tags=["panoramas"])
api_router.include_router(roi_router.router, prefix="/rois", tags=["rois"])
api_router.include_router(roi_point_router.router, prefix="/roi-points", tags=["roi-points"])
api_router.include_router(acquisition_router.router, prefix="/acquisitions", tags=["acquisitions"])
api_router.include_router(channel_router.router, prefix="/channels", tags=["channels"])
api_router.include_router(user_router.router, prefix="/users", tags=["users"])
api_router.include_router(dataset_router.router, prefix="/datasets", tags=["datasets"])
api_router.include_router(share_router.router, prefix="/share", tags=["share"])
api_router.include_router(analysis_router.router, prefix="/analysis", tags=["analysis"])
