import logging
import os

from sqlalchemy import Column, ForeignKey, Integer, String, Text, UniqueConstraint
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import Session, relationship

from app.core.utils import (
    FileUploadStatus,
    autocreate_directory_property,
    remove_location_upon_delete,
)
from app.db.base import CreatedAtMixin, DirectoryModel, MetaMixin
from app.modules.channel.db import Channel

logger = logging.getLogger(__name__)

#: Format string for acquisition locations
ACQUISITION_LOCATION_FORMAT = 'acquisition_{id}'


@remove_location_upon_delete
class Acquisition(DirectoryModel, MetaMixin, CreatedAtMixin):
    '''An *acquisition* contains all files belonging to one microscope image
    acquisition process. Note that in contrast to a *cycle*, an *acquisition*
    may contain more than one time point.

    The incentive to group files this way relates to the fact that most
    microscopes generate separate metadata files per *acquisition*.

    Attributes
    ----------
    channels: List[app.db_models.channel.Channel]
        channel files belonging to the acquisition
    '''

    __tablename__ = 'acquisition'
    __table_args__ = (UniqueConstraint('name', 'slide_id'),)

    #: str: name given by the user
    name = Column(String, index=True)

    #: int: acquisition width in pixels
    width = Column(Integer)

    #: int: acquisition height in pixels
    height = Column(Integer)

    #: str: description provided by the user
    description = Column(Text)

    #: int: ID of the parent slide
    slide_id = Column(
        Integer,
        ForeignKey('slide.id', onupdate='CASCADE', ondelete='CASCADE'),
        index=True
    )

    slide = relationship("Slide", back_populates="acquisitions")

    channels = relationship("Channel", back_populates="acquisition")

    def __init__(self, slide_id: int, name: str, width: int, height: int, description: str = '', meta: dict = None):
        '''
        Parameters
        ----------
        slide_id: int
            ID of the parent :class:`Slide <app.db_models.slide.Slide>`
        name: str
            name of the acquisition
        width: int
            width of the acquisition in pixels
        height: int
            height of the acquisition
        description: str, optional
            description of the acquisition
        meta: dict, optional
            meta data of the acquisition
        '''
        # TODO: ensure that name is unique within slide
        self.slide_id = slide_id
        self.name = name
        self.width = width
        self.height = height
        self.description = description
        self.meta = meta

    @hybrid_property
    def location(self) -> str:
        '''str: location were the acquisition content is stored'''
        if self._location is None:
            if self.id is None:
                raise AttributeError(
                    f'Acquisition "{self.name}" doesn\'t have an entry in the database yet. '
                    'Therefore, its location cannot be determined.'
                )
            self._location = os.path.join(
                self.slide.acquisitions_location,
                ACQUISITION_LOCATION_FORMAT.format(id=self.id)
            )
            if not os.path.exists(self._location):
                logger.debug(
                    'create location for acquisition "%s": %s',
                    self.name, self._location
                )
                os.mkdir(self._location)
        return self._location

    @location.setter
    def location(self, path_to_files: str):
        self._location = path_to_files

    @autocreate_directory_property
    def channels_location(self):
        '''str: location where channels files are stored'''
        return os.path.join(self.location, 'channels')

    @property
    def status(self) -> FileUploadStatus:
        '''str: upload status based on the status of channel files'''
        session = Session.object_session(self)
        channels = session.query(Channel.status). \
            filter_by(acquisition_id=self.id). \
            group_by(Channel.status). \
            all()
        if FileUploadStatus.UPLOADING in channels:
            return FileUploadStatus.UPLOADING
        elif len(channels) == 1 and FileUploadStatus.COMPLETE in channels:
            return FileUploadStatus.COMPLETE
        else:
            return FileUploadStatus.WAITING

    def to_dict(self) -> dict:
        '''Returns attributes as key-value pairs.

        Returns
        -------
        dict
        '''
        return {
            'id': self.id,
            'name': self.name,
            'description': self.description,
            'status': self.status
        }

    def json(self):
        return {
            'id': self.id,
            'slide_id': self.slide_id,
            'name': self.name,
            'description': self.description,
            'location': self.location,
            'meta': self.meta,
            'created_at': self.created_at,
            'channels': self.channels
        }

    def __repr__(self):
        return f'<Acquisition(id={self.id}, name={self.name})>'
