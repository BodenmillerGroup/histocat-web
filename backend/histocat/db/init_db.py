# make sure all SQL Alchemy models are imported before initializing DB
# otherwise, SQL Alchemy might fail to initialize properly relationships
# for more details: https://github.com/tiangolo/full-stack-fastapi-postgresql/issues/28

from histocat.config import config
from histocat.db.base import Base  # noqa
from histocat.modules.acquisition.models import AcquisitionModel  # noqa
from histocat.modules.dataset.models import DatasetModel  # noqa
from histocat.modules.experiment.models import ExperimentModel  # noqa
from histocat.modules.panorama.models import PanoramaModel  # noqa
from histocat.modules.share.models import ShareModel  # noqa
from histocat.modules.slide.models import SlideModel  # noqa
from histocat.modules.user import service
from histocat.modules.user.dto import UserCreateDto
from histocat.modules.user.models import UserModel  # noqa


def init_db(db_session):
    # Tables should be created with Alembic migrations
    # But if you don't want to use migrations, create
    # the tables un-commenting the next line
    # Base.metadata.create_all(bind=engine)

    user = service.get_by_email(db_session, email=config.FIRST_SUPERUSER)
    if not user:
        user_in = UserCreateDto(email=config.FIRST_SUPERUSER, password=config.FIRST_SUPERUSER_PASSWORD, is_admin=True, )
        user = service.create(db_session, params=user_in)
