import logging
import os
from shutil import rmtree

from cached_property import cached_property
from sqlalchemy import Column, Integer, ForeignKey, String, UniqueConstraint
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship, backref

from app.core.errors import DataError
from app.core.utils import remove_location_upon_delete, create_directory, autocreate_directory_property
from app.db.base import DirectoryModel, MetaMixin, CreatedAtMixin

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

    #: int: number of bytes used to encode intensity
    bit_depth = Column(Integer)

    #: int: maximum intensity value at which images get clipped at original
    #: bit depth before rescaling to 8-bit
    max_intensity = Column(Integer)

    #: int: minimum intensity value at which images get clipped at original
    #: bit depth before rescaling to 8-bit
    min_intensity = Column(Integer)

    #: int: ID of the parent acquisition
    acquisition_id = Column(
        Integer,
        ForeignKey('acquisition.id', onupdate='CASCADE', ondelete='CASCADE'),
        index=True
    )

    #: tmlib.models.experiment.Experiment: parent experiment
    acquisition = relationship(
        'Acquisition',
        backref=backref('channels', cascade='all, delete-orphan')
    )

    def __init__(self, name: str, metal: str, acquisition_id: int, bit_depth: int, meta: dict = None):
        '''
        Parameters
        ----------
        name: str
            name of the channel
        metal: str
            name of the metal
        acquisition_id: int
            ID of the parent
            :class:`Experiment <tmlib.models.experiment.Experiment>`
        bit_depth: int
            number of bits used to indicate intensity of pixels
        meta: dict, optional
            meta data of the channel
        '''
        self.name = name
        self.bit_depth = bit_depth
        self.acquisition_id = acquisition_id
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

    @cached_property
    def width(self) -> int:
        '''int: number of pixels along horizontal axis at highest resolution level
        '''
        logger.debug('retrieve channel "width" from parent acquisition')
        width = self.channel.acquisition.width
        if width is None:
            raise DataError('Channel width has not yet been calculated.')
        return width

    @cached_property
    def height(self) -> int:
        '''int: number of pixels along vertical axis at highest resolution level
        '''
        logger.debug('retrieve channel "height" from parent acquisition')
        height = self.channel.acquisition.height
        if height is None:
            raise DataError('Channel height has not yet been calculated.')
        return height

    def __repr__(self):
        return f'<Channel(id={self.id}, name={self.name})>'
