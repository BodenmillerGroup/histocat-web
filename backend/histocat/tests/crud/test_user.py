from fastapi.encoders import jsonable_encoder

from histocat.db.session import db_session
from histocat.modules.user import service
from histocat.modules.user.dto import UserCreateDto
from histocat.tests.utils.utils import random_lower_string


def test_create_user():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateDto(email=email, password=password)
    user = service.create(db_session, params=user_in)
    assert user.email == email
    assert hasattr(user, "password")


def test_authenticate_user():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateDto(email=email, password=password)
    user = service.create(db_session, params=user_in)
    authenticated_user = service.authenticate(db_session, email=email, password=password)
    assert authenticated_user
    assert user.email == authenticated_user.email


def test_not_authenticate_user():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user = service.authenticate(db_session, email=email, password=password)
    assert user is None


def test_check_if_user_is_active():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateDto(email=email, password=password)
    user = service.create(db_session, params=user_in)
    assert user.is_active is True


def test_check_if_user_is_active_inactive():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateDto(email=email, password=password, disabled=True)
    print(user_in)
    user = service.create(db_session, params=user_in)
    print(user)
    assert user.is_active


def test_check_if_user_is_admin():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateDto(email=email, password=password, is_admin=True)
    user = service.create(db_session, params=user_in)
    assert user.is_admin is True


def test_check_if_user_is_normal_user():
    username = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateDto(email=username, password=password)
    user = service.create(db_session, params=user_in)
    assert user.is_admin is False


def test_get_user():
    password = random_lower_string()
    username = f"{random_lower_string()}@test.com"
    user_in = UserCreateDto(email=username, password=password, is_admin=True)
    user = service.create(db_session, params=user_in)
    user_2 = service.get_by_id(db_session, user.id)
    assert user.email == user_2.email
    assert jsonable_encoder(user) == jsonable_encoder(user_2)
