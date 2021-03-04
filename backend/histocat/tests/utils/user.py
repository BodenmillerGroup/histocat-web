import requests

from histocat.config import config
from histocat.core.session import db_session
from histocat.core.user import service
from histocat.core.user.dto import UserCreateDto
from histocat.tests.utils.utils import random_lower_string


def user_authentication_headers(server_api, email, password):
    data = {"username": email, "password": password}

    r = requests.post(f"{server_api}{config.API_V1_STR}/auth/login", data=data)
    response = r.json()
    auth_token = response["access_token"]
    headers = {"Authorization": f"Bearer {auth_token}"}
    return headers


def create_random_user():
    email = random_lower_string()
    password = random_lower_string()
    user_in = UserCreateDto(username=email, email=email, password=password)
    user = service.create(db_session, params=user_in)
    return user
