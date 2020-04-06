from pydantic import BaseModel


class TokenDto(BaseModel):
    access_token: str
    token_type: str


class TokenPayloadDto(BaseModel):
    user_id: int
