import os

from pydantic import BaseSettings


class Settings(BaseSettings):
    API_V1_STR: str = "/api/v1"

    # absolute path to the directory where uploaded files are stored before import
    INBOX_DIRECTORY: str = "/data/inbox/"

    # absolute path to the root directory where all group data are located
    ROOT_DATA_DIRECTORY: str = "/data/groups/"

    SECRET_KEY: bytes = os.urandom(32)

    ACCESS_TOKEN_EXPIRE_MINUTES: int = 60 * 24 * 8  # 60 minutes * 24 hours * 8 days = 8 days

    SERVER_NAME: str
    SERVER_HOST: str

    # String of origins separated by commas, e.g: "http://localhost, http://localhost:4200, http://localhost:3000, http://localhost:8080, http://local.dockertoolbox.tiangolo.com"
    BACKEND_CORS_ORIGINS: str
    PROJECT_NAME: str

    POSTGRES_SERVER: str
    POSTGRES_USER: str
    POSTGRES_PASSWORD: str
    POSTGRES_DB: str

    SMTP_TLS: bool = True
    SMTP_PORT: int

    SMTP_HOST: str
    SMTP_USER: str
    SMTP_PASSWORD: str
    EMAILS_FROM_EMAIL: str
    EMAIL_RESET_TOKEN_EXPIRE_HOURS: int = 24
    EMAIL_CONFIRM_SIGNUP_EXPIRE_HOURS: int = 48
    EMAIL_TEMPLATES_DIR: str = "/app/histocat/email-templates"
    EMAILS_ENABLED: bool = True

    FIRST_SUPERUSER: str
    FIRST_SUPERUSER_PASSWORD: str

    USERS_OPEN_REGISTRATION: bool = False

    @property
    def SQLALCHEMY_DATABASE_URI(self):
        return f"postgresql://{self.POSTGRES_USER}:{self.POSTGRES_PASSWORD}@{self.POSTGRES_SERVER}/{self.POSTGRES_DB}"

    @property
    def EMAILS_FROM_NAME(self):
        return self.PROJECT_NAME


config = Settings()
