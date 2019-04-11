import logging
import os

from sqlalchemy import Column, String, Integer, Text, ForeignKey
from sqlalchemy import UniqueConstraint
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship, backref

from app.core.file_upload_status import FileUploadStatus
from app.core.utils import remove_location_upon_delete, autocreate_directory_property
from app.db_models.base import CreatedAtMixin, DirectoryModel, MetaMixin

logger = logging.getLogger(__name__)

#: Format string for slide locations
SLIDE_LOCATION_FORMAT = 'slide_{id}'


@remove_location_upon_delete
class Slide(DirectoryModel, MetaMixin, CreatedAtMixin):
    '''A *slide* represents a container with reservoirs for biological
    samples (referred to as *wells*).
    It's assumed that all images of a *slide* were acquired with the
    same microscope settings implying that each acquisition has the
    same number of *z-planes* and *channels*.

    Note
    ----
    For consistency, a *slide* is considered a single-well *plate*, i.e. a
    *plate* with only one *well*.

    Attributes
    ----------
    acquisitions: List[app.db_model.acquisition.Acqusition]
        acquisitions belonging to the slide
    '''

    __tablename__ = 'slide'
    __table_args__ = (UniqueConstraint('name'),)

    #: str: name given by user
    name = Column(String, index=True)

    #: str: description provided by user
    description = Column(Text)

    #: int: ID of parent experiment
    experiment_id = Column(
        Integer,
        ForeignKey('experiment.id', onupdate='CASCADE', ondelete='CASCADE'),
        index=True
    )

    #: tmlib.models.experiment.Experiment: parent experiment
    experiment = relationship(
        'Experiment',
        backref=backref('slides', cascade='all, delete-orphan')
    )

    def __init__(self, name: str, experiment_id: int, description: str = ''):
        '''
        Parameters
        ----------
        name: str
            name of the slide
        experiment_id: int
            ID of the parent
            :class:`Experiment <tmlib.models.experiment.Experiment>`
        description: str, optional
            description of the slide
        '''
        self.name = name
        self.description = description
        self.experiment_id = experiment_id

    @hybrid_property
    def location(self) -> str:
        '''str: location were the slide is stored'''
        if self._location is None:
            if self.id is None:
                raise AttributeError(
                    f'Slide "{self.name}" doesn\'t have an entry in the database yet. '
                    'Therefore, its location cannot be determined.'
                )
            self._location = os.path.join(
                self.experiment.slides_location,
                SLIDE_LOCATION_FORMAT.format(id=self.id)
            )
            if not os.path.exists(self._location):
                logger.debug(
                    'create location for slide "%s": %s',
                    self.name, self._location
                )
                os.mkdir(self._location)
        return self._location

    @location.setter
    def location(self, path_to_files: str):
        self._location = path_to_files

    @autocreate_directory_property
    def acquisitions_location(self) -> str:
        '''str: location where acquisitions are stored'''
        return os.path.join(self.location, 'acquisitions')

    @property
    def status(self) -> FileUploadStatus:
        '''str: upload status based on the status of acquisitions'''
        child_status = set([f.status for f in self.acquisitions])
        if FileUploadStatus.UPLOADING in child_status:
            return FileUploadStatus.UPLOADING
        elif len(child_status) == 1 and FileUploadStatus.COMPLETE in child_status:
            return FileUploadStatus.COMPLETE
        else:
            return FileUploadStatus.WAITING

    def __repr__(self):
        return f'<Slide(id={self.id}, name={self.name})>'
