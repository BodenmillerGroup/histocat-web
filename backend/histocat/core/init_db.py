from histocat.config import config
from histocat.core.acquisition.models import AcquisitionModel  # noqa
from histocat.core.base import Base  # noqa
from histocat.core.dataset.models import DatasetModel  # noqa
from histocat.core.gate.models import GateModel  # noqa
from histocat.core.group.models import GroupModel  # noqa
from histocat.core.member.models import MemberModel  # noqa
from histocat.core.model.models import ModelModel  # noqa
from histocat.core.panorama.models import PanoramaModel  # noqa
from histocat.core.pipeline.models import PipelineModel  # noqa
from histocat.core.preset.models import PresetModel  # noqa
from histocat.core.project.models import ProjectModel  # noqa
from histocat.core.result.models import ResultModel  # noqa
from histocat.core.security import get_password_hash
from histocat.core.slide.models import SlideModel  # noqa
from histocat.core.user import service
from histocat.core.user.dto import UserCreateDto
from histocat.core.user.models import UserModel  # noqa

# make sure all SQL Alchemy models are imported (app.db.base) before initializing DB
# otherwise, SQL Alchemy might fail to initialize relationships properly
# for more details: https://github.com/tiangolo/full-stack-fastapi-postgresql/issues/28


def init_db(db_session):
    # Tables should be created with Alembic migrations
    # But if you don't want to use migrations, create
    # the tables un-commenting the next line
    # Base.metadata.create_all(bind=engine)

    user = service.get_by_email(db_session, email=config.FIRST_SUPERUSER)
    if not user:
        hashed_password = get_password_hash(config.FIRST_SUPERUSER_PASSWORD)
        user_in = UserCreateDto(
            email=config.FIRST_SUPERUSER,
            password=hashed_password,
            is_admin=True,
        )
        user = service.create(db_session, params=user_in)
