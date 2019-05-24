import logging
import os
from shutil import rmtree

from sqlalchemy import Column, ForeignKey, Integer, String, UniqueConstraint
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship

from app.core.errors import DataError
from app.core.utils import (
    autocreate_directory_property,
    create_directory,
    remove_location_upon_delete,
)
from app.db.base import CreatedAtMixin, DirectoryModel, MetaMixin

logger = logging.getLogger(__name__)

#: Format string for channel locations
CHANNEL_LOCATION_FORMAT = 'channel_{id}'


@remove_location_upon_delete
class Channel(DirectoryModel, MetaMixin, CreatedAtMixin):
    '''A *channel* represents all *images* across different time points and
    spatial positions that were acquired with the same illumination and
    microscope filter settings.

    Attributes
    ----------
    image_files: List[tmlib.models.file.ChannelImageFile]
        images belonging to the channel
    '''

    __tablename__ = 'channel'
    __table_args__ = (UniqueConstraint('name'),)

    #: str: name given by the microscope or user
    name = Column(String, index=True)

    #: str: metal name
    metal = Column(String, index=True)

    #: int: metal mass
    mass = Column(Integer)

    #: int: maximum intensity value at which images get clipped at original bit depth before rescaling to 8-bit
    max_intensity = Column(Integer)

    #: int: minimum intensity value at which images get clipped at original bit depth before rescaling to 8-bit
    min_intensity = Column(Integer)

    #: int: ID of the parent acquisition
    acquisition_id = Column(
        Integer,
        ForeignKey('acquisition.id', onupdate='CASCADE', ondelete='CASCADE'),
        index=True
    )

    acquisition = relationship("Acquisition", back_populates="channels")

    def __init__(self, acquisition_id: int, name: str, metal: str, mass: int, max_intensity: int, min_intensity: int, meta: dict = None):
        '''
        Parameters
        ----------
        acquisition_id: int
            ID of the parent
            :class:`Experiment <tmlib.models.experiment.Experiment>`
        name: str
            name of the channel
        metal: str
            name of the metal
        mass: int
            mass of the metal
        max_intensity: int
            maximum intensity value
        min_intensity: int
            minimum intensity value
        meta: dict, optional
            meta data of the channel
        '''
        self.acquisition_id = acquisition_id
        self.name = name
        self.metal = metal
        self.mass = mass
        self.max_intensity = max_intensity
        self.min_intensity = min_intensity
        self.meta = meta

    @hybrid_property
    def location(self) -> str:
        '''str: location were channel content is stored'''
        if self._location is None:
            if self.id is None:
                raise AttributeError(
                    f'Channel "{self.name}" doesn\'t have an entry in the database yet. '
                    'Therefore, its location cannot be determined.'
                )
            self._location = os.path.join(
                self.acquisition.channels_location,
                CHANNEL_LOCATION_FORMAT.format(id=self.id)
            )
            if not os.path.exists(self._location):
                logger.debug(
                    'create location for channel "%s": %s',
                    self.name, self._location
                )
                create_directory(self._location)
        return self._location

    @autocreate_directory_property
    def image_files_location(self) -> str:
        '''str: location where image files are stored'''
        return os.path.join(self.location, 'images')

    def get_image_file_location(self, image_file_id: int) -> str:
        # TODO: It's not ideal to store them all in one directory. While modern
        # filesystems are able to handle this relatively well we should get
        # better performance using subdirectories.
        # Use a hash function to map image ID to subdirectory.
        return self.image_files_location

    def remove_image_files(self):
        '''Removes all image files on disk'''
        # TODO: walk sudirectories and delete all files
        rmtree(self.image_files_location)

    @autocreate_directory_property
    def illumstats_location(self):
        '''str: location where illumination statistics files are stored'''
        return os.path.join(self.location, 'illumstats')

    def width(self) -> int:
        '''int: number of pixels along horizontal axis at highest resolution level
        '''
        logger.debug('retrieve channel "width" from parent acquisition')
        width = self.acquisition.width
        if width is None:
            raise DataError('Channel width has not yet been calculated.')
        return width

    def height(self) -> int:
        '''int: number of pixels along vertical axis at highest resolution level
        '''
        logger.debug('retrieve channel "height" from parent acquisition')
        height = self.acquisition.height
        if height is None:
            raise DataError('Channel height has not yet been calculated.')
        return height

    def __repr__(self):
        return f'<Channel(id={self.id}, name={self.name})>'
