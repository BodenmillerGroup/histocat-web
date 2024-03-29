import logging
import os
import time
from datetime import datetime, timedelta
from io import BytesIO
from pathlib import Path
from shutil import rmtree
from typing import Optional

import dramatiq
import jwt
import sqlalchemy
from jwt.exceptions import InvalidTokenError

from histocat.config import config

password_reset_jwt_subject = "PasswordReset"
confirm_signup_jwt_subject = "ConfirmSignup"

logger = logging.getLogger(__name__)


def create_directory(location):
    """Creates a directory on disk in a safe way.

    Parameters
    ----------
    location: str
        absolute path to the directory that should be created
    """
    try:
        os.makedirs(location)
    except OSError as err:
        if err.errno != 17:
            raise


class autocreate_directory_property(object):
    """Decorator class that acts like a property.
    The value represents a path to a directory on disk. The directory is
    automatically created when it doesn't exist. Once created, the value
    is cached, so that there is no reattempt to create the directory.

    Raises
    ------
    TypeError
        when the value of the property doesn't have type basestring
    ValueError
        when the value of the property is empty

    Examples
    --------
    .. code:: python

        from tmlib.utils import autocreate_directory_property

        class Foo(object):

            @autocreate_directory_property
            def my_new_directory(self):
                return '/tmp/blabla'

        foo = Foo()
        foo.my_new_directory
    """

    def __init__(self, func):
        self.__doc__ = func.__doc__
        self.func = func

    def __get__(self, obj, cls):
        if obj is None:
            return self
        value = obj.__dict__[self.func.__name__] = self.func(obj)
        if not isinstance(value, str):
            raise TypeError('Value of property "%s" must have type string: %s' % (self.func.__name__, value))
        if not value:
            raise ValueError('Value of property "%s" cannot be empty.' % self.func.__name__)
        if not os.path.exists(value):
            logger.debug("create directory: %s", value)
            create_directory(value)
        return value


def delete_location(path: str):
    """Deletes a location on disk.

    Parameters
    ----------
    path: str
        absolute path to directory or file
    """
    if os.path.exists(path):
        logger.debug("remove location: %s", path)
        if os.path.isdir(path):
            rmtree(path)
        elif os.path.isfile(path):
            os.remove(path)


def remove_location_upon_delete(cls):
    """Decorator function for an database model class that
    automatically removes the `location` that represents an instance of the
    class on the filesystem once the corresponding row is deleted from the
    database table.

    Parameters
    ----------
    cls: tmlib.models.base.DeclarativeABCMeta
       implemenation of :class:`tmlib.models.base.FileSystemModel`

    Raises
    ------
    AttributeError
        when decorated class doesn't have a "location" attribute
    """

    def after_delete_callback(mapper, connection, target):
        delete_location(target.location)

    sqlalchemy.event.listen(cls, "after_delete", after_delete_callback)
    return cls


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if "log_time" in kw:
            name = kw.get("log_name", method.__name__.upper())
            kw["log_time"][name] = int((te - ts))
        else:
            print("%r  %2.2f s" % (method.__name__, (te - ts)))
        return result

    return timed


def send_test_email(email_to: str):
    project_name = config.PROJECT_NAME
    subject = f"{project_name} - Test email"
    with open(Path(config.EMAIL_TEMPLATES_DIR) / "test_email.html") as f:
        template_str = f.read()

    broker = dramatiq.get_broker()
    message = dramatiq.Message(
        actor_name="send_email",
        queue_name="email",
        args=(),
        kwargs={
            "email_to": email_to,
            "subject_template": subject,
            "html_template": template_str,
            "environment": {"project_name": config.PROJECT_NAME, "email": email_to},
        },
        options={},
    )
    broker.enqueue(message)


def send_new_account_email(email_to: str, username: str, password: str):
    project_name = config.PROJECT_NAME
    subject = f"{project_name} - New account for user {username}"
    with open(Path(config.EMAIL_TEMPLATES_DIR) / "new_account.html") as f:
        template_str = f.read()
    link = f"{config.PROTOCOL}://{config.DOMAIN}"

    broker = dramatiq.get_broker()
    message = dramatiq.Message(
        actor_name="send_email",
        queue_name="email",
        args=(),
        kwargs={
            "email_to": email_to,
            "subject_template": subject,
            "html_template": template_str,
            "environment": {
                "project_name": config.PROJECT_NAME,
                "username": username,
                "password": password,
                "email": email_to,
                "link": link,
            },
        },
        options={},
    )
    broker.enqueue(message)


def send_reset_password_email(email_to: str, email: str, token):
    project_name = config.PROJECT_NAME
    subject = f"{project_name} - Password recovery for user {email}"
    with open(Path(config.EMAIL_TEMPLATES_DIR) / "reset_password.html") as f:
        template_str = f.read()
    if hasattr(token, "decode"):
        use_token = token.decode()
    else:
        use_token = token
    link = f"{config.PROTOCOL}://{config.DOMAIN}/reset-password?token={use_token}"

    broker = dramatiq.get_broker()
    message = dramatiq.Message(
        actor_name="send_email",
        queue_name="email",
        args=(),
        kwargs={
            "email_to": email_to,
            "subject_template": subject,
            "html_template": template_str,
            "environment": {
                "project_name": config.PROJECT_NAME,
                "username": email,
                "email": email_to,
                "valid_hours": config.EMAIL_RESET_TOKEN_EXPIRE_HOURS,
                "link": link,
            },
        },
        options={},
    )
    broker.enqueue(message)


def send_confirm_signup_email(email_to: str, username: str, token):
    project_name = config.PROJECT_NAME
    subject = f"{project_name} - New account for user {username}"
    with open(Path(config.EMAIL_TEMPLATES_DIR) / "confirm_signup.html") as f:
        template_str = f.read()
    if hasattr(token, "decode"):
        use_token = token.decode()
    else:
        use_token = token
    link = f"{config.PROTOCOL}://{config.DOMAIN}/api/v1/auth/confirm-signup/{use_token}"

    broker = dramatiq.get_broker()
    message = dramatiq.Message(
        actor_name="send_email",
        queue_name="email",
        args=(),
        kwargs={
            "email_to": email_to,
            "subject_template": subject,
            "html_template": template_str,
            "environment": {
                "project_name": config.PROJECT_NAME,
                "username": username,
                "email": email_to,
                "valid_hours": config.EMAIL_CONFIRM_SIGNUP_EXPIRE_HOURS,
                "link": link,
            },
        },
        options={},
    )
    broker.enqueue(message)


def generate_password_reset_token(email):
    delta = timedelta(hours=config.EMAIL_RESET_TOKEN_EXPIRE_HOURS)
    now = datetime.utcnow()
    expires = now + delta
    exp = expires.timestamp()
    encoded_jwt = jwt.encode(
        {"exp": exp, "nbf": now, "sub": password_reset_jwt_subject, "email": email},
        config.JWT_SECRET,
        algorithm="HS256",
    )
    return encoded_jwt


def verify_password_reset_token(token) -> Optional[str]:
    try:
        decoded_token = jwt.decode(token, config.JWT_SECRET, algorithms=["HS256"])
        assert decoded_token["sub"] == password_reset_jwt_subject
        return decoded_token["email"]
    except InvalidTokenError:
        return None


def generate_confirm_signup_token(email: str):
    delta = timedelta(hours=config.EMAIL_CONFIRM_SIGNUP_EXPIRE_HOURS)
    now = datetime.utcnow()
    expires = now + delta
    exp = expires.timestamp()
    encoded_jwt = jwt.encode(
        {"exp": exp, "nbf": now, "sub": confirm_signup_jwt_subject, "email": email},
        config.JWT_SECRET,
        algorithm="HS256",
    )
    return encoded_jwt


def verify_confirm_signup_token(token) -> Optional[str]:
    try:
        decoded_token = jwt.decode(token, config.JWT_SECRET, algorithms=["HS256"])
        assert decoded_token["sub"] == confirm_signup_jwt_subject
        return decoded_token["email"]
    except InvalidTokenError:
        return None


async def stream_bytes(record: bytes, chunk_size: int = 65536):
    with BytesIO(record) as stream:
        data = stream.read(chunk_size)
        while data:
            yield data
            data = stream.read(chunk_size)
