import asyncio

import aioredis
import orjson
import redis

from histocat.config import config
from histocat.core.notifier import Message, notifier

UPDATES_CHANNEL_NAME = "updates"


class RedisManager:
    def __init__(self):
        self._pub = redis.Redis(host=config.REDIS_HOST, port=config.REDIS_PORT)
        self._sub: aioredis.Redis = None
        self._cache: aioredis.Redis = None

    @property
    def pub(self) -> redis.Redis:
        return self._pub

    @property
    def sub(self) -> aioredis.Redis:
        return self._sub

    @property
    def cache(self) -> aioredis.Redis:
        return self._cache

    async def start(self):
        self._pub = await aioredis.create_redis(f"redis://{config.REDIS_HOST}:{config.REDIS_PORT}")
        self._sub = await aioredis.create_redis(f"redis://{config.REDIS_HOST}:{config.REDIS_PORT}")
        self._cache = await aioredis.create_redis(f"redis://{config.REDIS_HOST}:{config.REDIS_PORT}")
        channels = await self.sub.subscribe(UPDATES_CHANNEL_NAME)
        updates_channel: aioredis.Channel = channels[0]
        asyncio.ensure_future(self._reader(updates_channel))

    async def stop(self):
        await self.sub.unsubscribe(UPDATES_CHANNEL_NAME)
        await self._close(self.pub)
        await self._close(self.sub)
        await self._close(self.cache)

    def publish(self, channel_name: str, message: Message):
        if self.pub is not None:
            self.pub.publish(channel_name, orjson.dumps(message.to_json()))

    # async def publish_async(self, channel_name: str, message: Message):
    #     if self.pub is not None:
    #         await self.pub.publish_json(channel_name, message.to_json())

    async def _reader(self, channel: aioredis.Channel):
        while await channel.wait_message():
            json = await channel.get_json()
            message = Message.from_json(json)
            await notifier.push(message)

    async def _close(self, redis: aioredis.Redis):
        if redis is not None:
            redis.close()
            await redis.wait_closed()


redis_manager = RedisManager()
