from redis import Redis

from app.core.notifier import Message

redis = Redis(host="redis")
pubsub = redis.pubsub()


def publish(channel: str, message: Message):
    redis.publish(channel, message)


def subscribe(channel: str, handler):
    pubsub.subscribe(**{channel: handler})
