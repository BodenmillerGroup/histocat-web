from fastapi.encoders import jsonable_encoder

from app.db.session import db_session
from app.modules.user.models import UserInCreateModel
from app.tests.utils.utils import random_lower_string
from app.modules.user import crud


def test_create_user():
    email = random_lower_string()
    password = random_lower_string()
    user_in = UserInCreateModel(email=email, password=password)
    user = crud.create(db_session, params=user_in)
    assert user.email == email
    assert hasattr(user, "hashed_password")


def test_authenticate_user():
    email = random_lower_string()
    password = random_lower_string()
    user_in = UserInCreateModel(email=email, password=password)
    user = crud.create(db_session, params=user_in)
    authenticated_user = crud.authenticate(
        db_session, email=email, password=password
    )
    assert authenticated_user
    assert user.email == authenticated_user.email


def test_not_authenticate_user():
    email = random_lower_string()
    password = random_lower_string()
    user = crud.authenticate(db_session, email=email, password=password)
    assert user is None


def test_check_if_user_is_active():
    email = random_lower_string()
    password = random_lower_string()
    user_in = UserInCreateModel(email=email, password=password)
    user = crud.create(db_session, params=user_in)
    is_active = crud.is_active(user)
    assert is_active is True


def test_check_if_user_is_active_inactive():
    email = random_lower_string()
    password = random_lower_string()
    user_in = UserInCreateModel(email=email, password=password, disabled=True)
    print(user_in)
    user = crud.create(db_session, params=user_in)
    print(user)
    is_active = crud.is_active(user)
    print(is_active)
    assert is_active


def test_check_if_user_is_superuser():
    email = random_lower_string()
    password = random_lower_string()
    user_in = UserInCreateModel(email=email, password=password, is_superuser=True)
    user = crud.create(db_session, params=user_in)
    is_superuser = crud.is_superuser(user)
    assert is_superuser is True


def test_check_if_user_is_superuser_normal_user():
    username = random_lower_string()
    password = random_lower_string()
    user_in = UserInCreateModel(email=username, password=password)
    user = crud.create(db_session, params=user_in)
    is_superuser = crud.is_superuser(user)
    assert is_superuser is False


def test_get_user():
    password = random_lower_string()
    username = random_lower_string()
    user_in = UserInCreateModel(email=username, password=password, is_superuser=True)
    user = crud.create(db_session, params=user_in)
    user_2 = crud.get(db_session, id=user.id)
    assert user.email == user_2.email
    assert jsonable_encoder(user) == jsonable_encoder(user_2)
