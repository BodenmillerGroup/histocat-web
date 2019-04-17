from fastapi import APIRouter

from app.modules.auth import router as auth_router
from app.modules.core import router as core_router
from app.modules.experiment import router as experiment_router
from app.modules.slide import router as slide_router
from app.modules.acquisition import router as acquisition_router
from app.modules.channel import router as channel_router
from app.modules.user import router as user_router

api_router = APIRouter()
api_router.include_router(auth_router.router)
api_router.include_router(core_router.router)
api_router.include_router(experiment_router.router)
api_router.include_router(slide_router.router)
api_router.include_router(acquisition_router.router)
api_router.include_router(channel_router.router)
api_router.include_router(user_router.router)
