# make sure all SQL Alchemy models are imported before initializing DB
# otherwise, SQL Alchemy might fail to initialize properly relationships
# for more details: https://github.com/tiangolo/full-stack-fastapi-postgresql/issues/28

from app.db.base import Base  # noqa
from app.modules.user.db import User  # noqa
from app.modules.experiment.db import Experiment  # noqa
from app.modules.slide.db import Slide  # noqa
from app.modules.panorama.db import Panorama  # noqa
from app.modules.roi.db import ROI  # noqa
from app.modules.roi_point.db import ROIPoint  # noqa
from app.modules.acquisition.db import Acquisition  # noqa
from app.modules.channel.db import Channel  # noqa
from app.modules.dataset.db import Dataset  # noqa
from app.modules.share.db import Share  # noqa
from app.modules.acquisition_artifact.db import AcquisitionArtifact  # noqa

from app.core import config
from app.modules.user import crud
from app.modules.user.models import UserCreateModel


def init_db(db_session):
    # Tables should be created with Alembic migrations
    # But if you don't want to use migrations, create
    # the tables un-commenting the next line
    # Base.metadata.create_all(bind=engine)

    user = crud.get_by_email(db_session, email=config.FIRST_SUPERUSER)
    if not user:
        user_in = UserCreateModel(
            email=config.FIRST_SUPERUSER,
            password=config.FIRST_SUPERUSER_PASSWORD,
            is_superuser=True,
        )
        user = crud.create(db_session, params=user_in)
