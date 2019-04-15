from pydantic import BaseModel


class MsgModel(BaseModel):
    msg: str
