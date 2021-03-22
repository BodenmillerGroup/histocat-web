from typing import Optional

from sqlalchemy.orm import Session
from sqlalchemy.orm.attributes import flag_modified

from histocat.core.project import service as project_service

from ..slide.models import SlideModel
from .dto import AcquisitionCreateDto, ChannelUpdateDto
from .models import AcquisitionModel


def get_by_id(session: Session, id: int) -> Optional[AcquisitionModel]:
    return session.query(AcquisitionModel).filter(AcquisitionModel.id == id).first()


def get_by_origin_id(session: Session, *, slide_id: int, origin_id: int) -> Optional[AcquisitionModel]:
    return (
        session.query(AcquisitionModel)
        .filter(AcquisitionModel.slide_id == slide_id, AcquisitionModel.origin_id == origin_id)
        .first()
    )


def get_by_project_id_and_description(
    session: Session, *, project_id: int, description: str
) -> Optional[AcquisitionModel]:
    return (
        session.query(AcquisitionModel)
        .join(SlideModel.acquisitions)
        .filter(SlideModel.project_id == project_id, AcquisitionModel.description == description)
        .first()
    )


def create(session: Session, params: AcquisitionCreateDto) -> AcquisitionModel:
    data = params.dict()
    item = AcquisitionModel(**data)
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def update_custom_label(session: Session, *, item: AcquisitionModel, params: ChannelUpdateDto):
    for a in item.slide.acquisitions:
        channel_data = a.channels.get(params.name)
        channel_data["customLabel"] = params.customLabel
        a.channels[params.name] = channel_data

        # TODO: https://stackoverflow.com/questions/42559434/updates-to-json-field-dont-persist-to-db
        flag_modified(a, "channels")
        session.add(a)
        session.commit()
        session.refresh(a)
    return project_service.get_data(session, id=item.slide.project_id)


def delete_by_id(session: Session, id: int):
    item = get_by_id(session, id)
    session.delete(item)
    session.commit()
    return item
