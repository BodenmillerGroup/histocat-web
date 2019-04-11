import logging
import os

from sqlalchemy import Column, String, Integer, Text, ForeignKey, UniqueConstraint
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship, backref, Session

from app.core.file_upload_status import FileUploadStatus
from app.core.utils import remove_location_upon_delete, autocreate_directory_property
from app.db_models.base import DirectoryModel, MetaMixin, CreatedAtMixin
from app.db_models.channel import Channel

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

    #: str: description provided by the user
    description = Column(Text)

    #: int: number of pixels along *x*-axis
    width = Column(Integer)

    #: int: number of pixels along *y*-axis
    height = Column(Integer)

    #: int: ID of the parent slide
    slide_id = Column(
        Integer,
        ForeignKey('slide.id', onupdate='CASCADE', ondelete='CASCADE'),
        index=True
    )

    #: app.db_models.slide.Slide: slide to which the acquisition belongs
    slide = relationship(
        'Slide',
        backref=backref('acquisitions', cascade='all, delete-orphan')
    )

    def __init__(self, name: str, slide_id: int, description: str = ''):
        '''
        Parameters
        ----------
        name: str
            name of the acquisition
        slide_id: int
            ID of the parent :class:`Slide <app.db_models.slide.Slide>`
        description: str, optional
            description of the acquisition
        '''
        # TODO: ensure that name is unique within slide
        self.name = name
        self.description = description
        self.slide_id = slide_id

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
                self.plate.acquisitions_location,
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

    def __repr__(self):
        return f'<Acquisition(id={self.id}, name={self.name})>'
