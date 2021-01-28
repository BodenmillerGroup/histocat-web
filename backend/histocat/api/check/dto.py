from pydantic import BaseModel


class MsgDto(BaseModel):
    msg: str
