import logging
import os
from datetime import datetime
from shutil import rmtree

from sqlalchemy import Column, ForeignKey, Integer, String, UniqueConstraint, DateTime
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship
from sqlalchemy.sql.functions import now

from app.core.errors import DataError
from app.core.utils import (
    autocreate_directory_property,
    remove_location_upon_delete,
)
from app.db.base import Base

logger = logging.getLogger(__name__)

#: Format string for channel locations
CHANNEL_LOCATION_FORMAT = "channel_{id}"


@remove_location_upon_delete
class Channel(Base):
    """
    A *channel* represents *image* that was acquired with the same illumination and microscope filter settings.
    """

    __tablename__ = "channel"
    __table_args__ = (UniqueConstraint("name"),)

    id: int = Column(Integer, primary_key=True, index=True)
    #: ID of the parent acquisition
    acquisition_id: int = Column(
        Integer,
        ForeignKey("acquisition.id", onupdate="CASCADE", ondelete="CASCADE"),
        index=True,
    )
    #: location of the channel, e.g. absolute path to a directory on disk
    location: str = Column("location", String(4096))
    #: name
    name: str = Column(String, index=True)
    #: metal name
    metal: str = Column(String, index=True)
    #: metal mass
    mass: int = Column(Integer)
    #: maximum intensity value at which images get clipped at original bit depth before rescaling to 8-bit
    max_intensity: int = Column(Integer)
    #: minimum intensity value at which images get clipped at original bit depth before rescaling to 8-bit
    min_intensity: int = Column(Integer)
    meta: dict = Column(JSONB)
    created_at: datetime = Column(DateTime, default=now(), nullable=False)

    acquisition = relationship("Acquisition", back_populates="channels")

    def __init__(
        self,
        acquisition_id: int,
        name: str,
        metal: str,
        mass: int,
        max_intensity: int,
        min_intensity: int,
        meta: dict = None,
    ):
        self.acquisition_id = acquisition_id
        self.name = name
        self.metal = metal
        self.mass = mass
        self.max_intensity = max_intensity
        self.min_intensity = min_intensity
        self.meta = meta

    @autocreate_directory_property
    def image_files_location(self) -> str:
        """
        Location where image files are stored
        """
        return os.path.join(self.location, "images")

    def get_image_file_location(self, image_file_id: int) -> str:
        # TODO: It's not ideal to store them all in one directory. While modern
        # filesystems are able to handle this relatively well we should get
        # better performance using subdirectories.
        # Use a hash function to map image ID to subdirectory.
        return self.image_files_location

    def remove_image_files(self):
        """
        Removes all image files on disk
        """
        # TODO: walk sudirectories and delete all files
        rmtree(self.image_files_location)

    @autocreate_directory_property
    def illumstats_location(self) -> str:
        """
        Location where illumination statistics files are stored
        """
        return os.path.join(self.location, "illumstats")

    def width(self) -> int:
        """
        Number of pixels along horizontal axis at highest resolution level
        """
        logger.debug('retrieve channel "width" from parent acquisition')
        width = self.acquisition.width
        if width is None:
            raise DataError("Channel width has not yet been calculated.")
        return width

    def height(self) -> int:
        """
        Number of pixels along vertical axis at highest resolution level
        """
        logger.debug('retrieve channel "height" from parent acquisition')
        height = self.acquisition.height
        if height is None:
            raise DataError("Channel height has not yet been calculated.")
        return height

    def __repr__(self):
        return f"<Channel(id={self.id}, name={self.name})>"
