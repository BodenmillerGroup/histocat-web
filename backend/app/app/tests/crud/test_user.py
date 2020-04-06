from fastapi.encoders import jsonable_encoder

from app.db.session import db_session
from app.modules.user import crud
from app.modules.user.models import UserCreateModel
from app.tests.utils.utils import random_lower_string


def test_create_user():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateModel(email=email, password=password)
    user = crud.create(db_session, params=user_in)
    assert user.email == email
    assert hasattr(user, "password")


def test_authenticate_user():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateModel(email=email, password=password)
    user = crud.create(db_session, params=user_in)
    authenticated_user = crud.authenticate(db_session, email=email, password=password)
    assert authenticated_user
    assert user.email == authenticated_user.email


def test_not_authenticate_user():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user = crud.authenticate(db_session, email=email, password=password)
    assert user is None


def test_check_if_user_is_active():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateModel(email=email, password=password)
    user = crud.create(db_session, params=user_in)
    assert user.is_active is True


def test_check_if_user_is_active_inactive():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateModel(email=email, password=password, disabled=True)
    print(user_in)
    user = crud.create(db_session, params=user_in)
    print(user)
    assert user.is_active


def test_check_if_user_is_admin():
    email = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateModel(email=email, password=password, is_admin=True)
    user = crud.create(db_session, params=user_in)
    assert user.is_admin is True


def test_check_if_user_is_normal_user():
    username = f"{random_lower_string()}@test.com"
    password = random_lower_string()
    user_in = UserCreateModel(email=username, password=password)
    user = crud.create(db_session, params=user_in)
    assert user.is_admin is False


def test_get_user():
    password = random_lower_string()
    username = f"{random_lower_string()}@test.com"
    user_in = UserCreateModel(email=username, password=password, is_admin=True)
    user = crud.create(db_session, params=user_in)
    user_2 = crud.get(db_session, id=user.id)
    assert user.email == user_2.email
    assert jsonable_encoder(user) == jsonable_encoder(user_2)
