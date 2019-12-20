from dataclasses import dataclass
from typing import Any, Dict, List

from starlette.websockets import WebSocket


@dataclass
class Message:
    experiment_id: int
    type: str
    payload: Any = None

    @staticmethod
    def from_json(json: Dict[str, Any]):
        return Message(json.get("experiment_id"), json.get("type"), json.get("payload"))

    def to_json(self):
        return {
            "experiment_id": self.experiment_id,
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

    async def connect(self, websocket: WebSocket, experiment_id: int):
        await websocket.accept()
        if experiment_id not in self.connections:
            self.connections[experiment_id] = []
        self.connections.get(experiment_id).append(websocket)

    def remove(self, websocket: WebSocket, experiment_id: int):
        if experiment_id in self.connections:
            self.connections.get(experiment_id).remove(websocket)

    async def _notify(self, message: Message):
        living_connections = []
        json = message.to_json()
        if message.experiment_id in self.connections:
            while len(self.connections.get(message.experiment_id)) > 0:
                # Looping like this is necessary in case a disconnection is handled
                # during await websocket.send_text(message)
                websocket = self.connections.get(message.experiment_id).pop()
                await websocket.send_json(json)
                living_connections.append(websocket)
            self.connections[message.experiment_id] = living_connections


notifier = Notifier()
