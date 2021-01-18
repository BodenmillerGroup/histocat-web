import logging
import os

import dramatiq
from dramatiq.brokers.rabbitmq import RabbitmqBroker
from fastapi import FastAPI
from starlette.middleware.cors import CORSMiddleware
from starlette.websockets import WebSocket, WebSocketDisconnect

from histocat.api.router import api_router
from histocat.config import config
from histocat.core.notifier import notifier
from histocat.core.redis_manager import redis_manager

rabbitmq_broker = RabbitmqBroker(host="rabbitmq", connection_attempts=10)
dramatiq.set_broker(rabbitmq_broker)

logger = logging.getLogger(__name__)


if os.environ.get("BACKEND_ENV") == "development":
    try:
        # VS Code Debugging
        # Allow other computers to attach to ptvsd at this IP address and port.
        import ptvsd

        ptvsd.enable_attach(address=("0.0.0.0", 5678), redirect_output=True)

        # PyCharm Debugging
        # TODO: Don't forget to modify IP address!!
        # import pydevd_pycharm
        # pydevd_pycharm.settrace('192.168.1.129', port=5679, stdoutToServer=True, stderrToServer=True, suspend=False)

    except Exception as e:
        logger.error(e)


app = FastAPI(
    title=config.PROJECT_NAME, openapi_url=f"{config.API_V1_STR}/openapi.json", docs_url="/docs", redoc_url="/redoc",
)

# GZIP compression
# app.add_middleware(GZipMiddleware, minimum_size=100000)

# CORS
origins = []

# Set all CORS enabled origins
if config.BACKEND_CORS_ORIGINS:
    origins_raw = config.BACKEND_CORS_ORIGINS.split(",")
    for origin in origins_raw:
        use_origin = origin.strip()
        origins.append(use_origin)
    app.add_middleware(
        CORSMiddleware, allow_origins=origins, allow_credentials=True, allow_methods=["*"], allow_headers=["*"],
    )

app.include_router(api_router, prefix=config.API_V1_STR)


@app.websocket("/ws/{project_id}")
async def project_websocket_endpoint(websocket: WebSocket, project_id: int, token: str = None):
    if not token:
        raise Exception("WebSocket authorization token is missing")

    await notifier.connect(websocket, project_id)
    try:
        while True:
            await websocket.receive_text()
    except WebSocketDisconnect:
        notifier.remove(websocket, project_id)


@app.on_event("startup")
async def startup():
    await notifier.start()
    await redis_manager.start()


@app.on_event("shutdown")
async def shutdown():
    await notifier.stop()
    await redis_manager.stop()
