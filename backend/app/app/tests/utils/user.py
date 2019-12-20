import requests

from app.core import config
from app.db.session import db_session
from app.modules.user import crud
from app.modules.user.models import UserCreateModel
from app.tests.utils.utils import random_lower_string


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
    user_in = UserCreateModel(username=email, email=email, password=password)
    user = crud.create(db_session, params=user_in)
    return user
