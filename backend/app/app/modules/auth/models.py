from pydantic import BaseModel


class TokenModel(BaseModel):
    access_token: str
    token_type: str


class TokenPayloadModel(BaseModel):
    user_id: int
