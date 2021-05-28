import uuid

from pydantic import BaseSettings


class Settings(BaseSettings):
    API_V1_STR: str = "/api/v1"

    # absolute path to the directory where uploaded files are stored before import
    INBOX_DIRECTORY: str = "/data/inbox/"

    # absolute path to the directory where uploaded models are stored
    MODELS_DIRECTORY: str = "/data/models/"

    # absolute path to the root directory where all group data are located
    ROOT_DATA_DIRECTORY: str = "/data/groups/"

    JWT_SECRET: str = str(uuid.uuid4())

    ACCESS_TOKEN_EXPIRE_MINUTES: int = 60 * 24 * 8  # 60 minutes * 24 hours * 8 days = 8 days

    PROTOCOL: str
    DOMAIN: str

    PROJECT_NAME: str

    POSTGRES_SERVER: str
    POSTGRES_USER: str
    POSTGRES_PASSWORD: str
    POSTGRES_DB: str

    REDIS_HOST: str
    REDIS_PORT: int

    SMTP_TLS: bool = True
    SMTP_HOST: str
    SMTP_PORT: int
    SMTP_USER: str
    SMTP_PASSWORD: str

    EMAILS_FROM: str
    EMAIL_RESET_TOKEN_EXPIRE_HOURS: int = 24
    EMAIL_CONFIRM_SIGNUP_EXPIRE_HOURS: int = 48
    EMAIL_TEMPLATES_DIR: str = "/app/histocat/api/email-templates"
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
