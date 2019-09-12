import logging
import os

from fastapi import FastAPI
from starlette.middleware.cors import CORSMiddleware
from starlette.requests import Request
from starlette.websockets import WebSocket, WebSocketDisconnect

from app.api.api_v1.api import api_router
from app.core import config
from app.core.notifier import notifier
from app.core.redis_manager import redis_manager
from app.db.session import Session

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


if os.environ.get("BACKEND_ENV") == "development":
    try:
        # VS Code Debugging
        # Allow other computers to attach to ptvsd at this IP address and port.
        # import ptvsd
        # ptvsd.enable_attach(address=('0.0.0.0', 5678), redirect_output=True)

        # PyCharm Debugging
        # TODO: Don't forget to modify IP address!!
        # import pydevd_pycharm
        # pydevd_pycharm.settrace('130.60.106.31', port=5679, stdoutToServer=True, stderrToServer=True, suspend=False)

        pass
    except Exception as e:
        logger.error(e)

app = FastAPI(
    title=config.PROJECT_NAME,
    openapi_url=f"{config.API_V1_STR}/openapi.json",
)

# CORS
origins = []

# Set all CORS enabled origins
if config.BACKEND_CORS_ORIGINS:
    origins_raw = config.BACKEND_CORS_ORIGINS.split(",")
    for origin in origins_raw:
        use_origin = origin.strip()
        origins.append(use_origin)
    app.add_middleware(
        CORSMiddleware,
        allow_origins=origins,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    ),

app.include_router(api_router, prefix=config.API_V1_STR)


@app.middleware("http")
async def db_session_middleware(request: Request, call_next):
    request.state.db = Session()
    response = await call_next(request)
    request.state.db.close()
    return response


@app.websocket("/ws/{experiment_id}")
async def experiment_websocket_endpoint(
    websocket: WebSocket,
    experiment_id: int,
    token: str = None,
):
    if not token:
        raise Exception("WebSocket authorization token is missing")

    await notifier.connect(websocket, experiment_id)
    try:
        while True:
            await websocket.receive_text()
    except WebSocketDisconnect:
        notifier.remove(websocket, experiment_id)


@app.on_event("startup")
async def startup():
    await notifier.start()
    await redis_manager.start()


@app.on_event("shutdown")
async def shutdown():
    await notifier.stop()
    await redis_manager.stop()
