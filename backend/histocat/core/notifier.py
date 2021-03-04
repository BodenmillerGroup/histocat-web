from dataclasses import dataclass
from typing import Any, Dict, List

from starlette.websockets import WebSocket


@dataclass
class Message:
    project_id: int
    type: str
    payload: Any = None

    @staticmethod
    def from_json(json: Dict[str, Any]):
        return Message(json.get("project_id"), json.get("type"), json.get("payload"))

    def to_json(self):
        return {
            "project_id": self.project_id,
            "type": self.type,
            "payload": self.payload,
        }


class Notifier:
    def __init__(self):
        self.connections: Dict[int, List[WebSocket]] = {}
        self.generator = self.get_notification_generator()

    async def start(self):
        # Prime the push notification generator
        await self.generator.asend(None)

    async def stop(self):
        await self.generator.aclose()

    async def get_notification_generator(self):
        while True:
            message = yield
            await self._notify(message)

    async def push(self, message: Message):
        await self.generator.asend(message)

    async def connect(self, websocket: WebSocket, project_id: int):
        await websocket.accept()
        if project_id not in self.connections:
            self.connections[project_id] = []
        self.connections.get(project_id).append(websocket)

    def remove(self, websocket: WebSocket, project_id: int):
        if project_id in self.connections:
            self.connections.get(project_id).remove(websocket)

    async def _notify(self, message: Message):
        living_connections = []
        json = message.to_json()
        if message.project_id in self.connections:
            while len(self.connections.get(message.project_id)) > 0:
                # Looping like this is necessary in case a disconnection is handled
                # during await websocket.send_text(message)
                websocket = self.connections.get(message.project_id).pop()
                await websocket.send_json(json)
                living_connections.append(websocket)
            self.connections[message.project_id] = living_connections


notifier = Notifier()
