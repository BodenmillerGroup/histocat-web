import requests

from histocat.core import config
from histocat.db.session import db_session
from histocat.modules.user import service
from histocat.modules.user.dto import UserCreateDto
from histocat.tests.utils.user import user_authentication_headers
from histocat.tests.utils.utils import get_server_api, random_lower_string


def test_get_users_admin_me(superuser_token_headers):
    server_api = get_server_api()
    r = requests.get(f"{server_api}{config.API_V1_STR}/users/profile", headers=superuser_token_headers)
    current_user = r.json()
    assert current_user
    assert current_user["is_active"] is True
    assert current_user["is_admin"]
    assert current_user["email"] == config.FIRST_SUPERUSER


def test_create_user_new_email(superuser_token_headers):
    server_api = get_server_api()
    username = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    data = {"email": username, "password": password}
    r = requests.post(f"{server_api}{config.API_V1_STR}/users/", headers=superuser_token_headers, json=data,)
    assert 200 <= r.status_code < 300
    created_user = r.json()
    user = service.get_by_email(db_session, email=username)
    assert user.email == created_user["email"]


def test_get_existing_user(superuser_token_headers):
    server_api = get_server_api()
    username = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateDto(email=username, password=password)
    user = service.create(db_session, params=user_in)
    user_id = user.id
    r = requests.get(f"{server_api}{config.API_V1_STR}/users/{user_id}", headers=superuser_token_headers,)
    assert 200 <= r.status_code < 300
    api_user = r.json()
    user = service.get_by_email(db_session, email=username)
    assert user.email == api_user["email"]


def test_create_user_existing_username(superuser_token_headers):
    server_api = get_server_api()
    username = f"{random_lower_string()}@test.com"
    # username = email
    password = random_lower_string()
    user_in = UserCreateDto(email=username, password=password)
    user = service.create(db_session, params=user_in)
    data = {"email": username, "password": password}
    r = requests.post(f"{server_api}{config.API_V1_STR}/users/", headers=superuser_token_headers, json=data,)
    created_user = r.json()
    assert r.status_code == 400
    assert "_id" not in created_user


def test_create_user_by_normal_user():
    server_api = get_server_api()
    username = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateDto(email=username, password=password)
    user = service.create(db_session, params=user_in)
    user_token_headers = user_authentication_headers(server_api, username, password)
    data = {"email": username, "password": password}
    r = requests.post(f"{server_api}{config.API_V1_STR}/users/", headers=user_token_headers, json=data)
    assert r.status_code == 400


def test_retrieve_users(superuser_token_headers):
    server_api = get_server_api()
    username = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateDto(email=username, password=password)
    user = service.create(db_session, params=user_in)

    username2 = f"{random_lower_string()}@test.com"
    password2 = random_lower_string()
    user_in2 = UserCreateDto(email=username2, password=password2)
    user2 = service.create(db_session, params=user_in2)

    r = requests.get(f"{server_api}{config.API_V1_STR}/users/", headers=superuser_token_headers)
    all_users = r.json()

    assert len(all_users) > 1
    for user in all_users:
        assert "email" in user
