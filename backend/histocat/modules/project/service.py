import logging
import os
from typing import List, Optional, Set

from sqlalchemy import func
from sqlalchemy.orm import Session

from .dto import ProjectCreateDto, ProjectUpdateDto
from .models import ProjectModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[ProjectModel]:
    return session.query(ProjectModel).filter(ProjectModel.id == id).first()


def get_by_name(session: Session, *, name: str) -> Optional[ProjectModel]:
    return session.query(ProjectModel).filter(ProjectModel.name == name).first()


def get_group_projects(session: Session, *, group_id: int) -> List[Optional[ProjectModel]]:
    return session.query(ProjectModel).order_by(ProjectModel.id.desc()).filter(ProjectModel.group_id == group_id).all()


def get_tags(session: Session, *, group_id: int) -> Set[str]:
    items = session.query(func.unnest(ProjectModel.tags)).filter(ProjectModel.group_id == group_id).distinct().all()
    return {e[0] for e in items}


def create(session: Session, *, group_id: int, params: ProjectCreateDto, member_id: int) -> ProjectModel:
    entity = ProjectModel(
        group_id=group_id, name=params.name, description=params.description, tags=params.tags, member_id=member_id
    )
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(entity.group.projects_location, str(entity.id))
    if not os.path.exists(entity.location):
        logger.debug(f"Create location for project {entity.id}: {entity.location}")
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)

    return entity


def update(session: Session, *, item: ProjectModel, params: ProjectUpdateDto) -> ProjectModel:
    data = item.as_dict()
    update_data = params.dict(exclude_unset=True)
    for field in data:
        if field in update_data:
            setattr(item, field, update_data[field])
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def remove(session: Session, *, id: int):
    item = session.query(ProjectModel).filter(ProjectModel.id == id).first()
    session.delete(item)
    session.commit()
    return item


def get_data(session: Session, *, id: int) -> Optional[ProjectModel]:
    result = session.query(ProjectModel).filter(ProjectModel.id == id).first()
    return result
