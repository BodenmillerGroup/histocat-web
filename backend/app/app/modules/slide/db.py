import logging
import os

from sqlalchemy import Column, ForeignKey, Integer, String, Text, UniqueConstraint
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship

from app.core.utils import (
    FileUploadStatus,
    autocreate_directory_property,
    remove_location_upon_delete,
)
from app.db.base import CreatedAtMixin, DirectoryModel, MetaMixin

logger = logging.getLogger(__name__)

#: Format string for slide locations
SLIDE_LOCATION_FORMAT = "slide_{id}"


@remove_location_upon_delete
class Slide(DirectoryModel, MetaMixin, CreatedAtMixin):
    """A *slide* represents a container with reservoirs for biological
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
    """

    __tablename__ = "slide"
    __table_args__ = (UniqueConstraint("name"),)

    #: str: name given by user
    name = Column(String, index=True)

    #: str: original slide filename
    filename = Column(String(4096))

    #: int: slide width in μm
    width_um = Column(Integer)

    #: int: slide height in μm
    height_um = Column(Integer)

    #: str: description provided by user
    description = Column(Text)

    #: int: ID of parent experiment
    experiment_id = Column(
        Integer,
        ForeignKey("experiment.id", onupdate="CASCADE", ondelete="CASCADE"),
        index=True,
    )

    experiment = relationship("Experiment", back_populates="slides")

    acquisitions = relationship("Acquisition", back_populates="slide")

    def __init__(
        self,
        experiment_id: int,
        name: str,
        filename: str,
        width_um: int,
        height_um: int,
        description: str = "",
        meta: dict = None,
    ):
        """
        Parameters
        ----------
        experiment_id: int
            ID of the parent
            :class:`Experiment <tmlib.models.experiment.Experiment>`
        name: str
            name of the slide
        filename: str
            name of the original slide file
        width_um: int
            slide width in μm
        height_um: int
            slide height in μm
        description: str, optional
            description of the slide
        meta: dict, optional
            meta data of the slide
        """
        self.experiment_id = experiment_id
        self.name = name
        self.filename = filename
        self.width_um = width_um
        self.height_um = height_um
        self.description = description
        self.meta = meta

    @hybrid_property
    def location(self) -> str:
        """str: location were the slide is stored"""
        if self._location is None:
            if self.id is None:
                raise AttributeError(
                    f'Slide "{self.name}" doesn\'t have an entry in the database yet. '
                    "Therefore, its location cannot be determined."
                )
            self._location = os.path.join(
                self.experiment.slides_location,
                SLIDE_LOCATION_FORMAT.format(id=self.id),
            )
            if not os.path.exists(self._location):
                logger.debug(
                    'create location for slide "%s": %s', self.name, self._location
                )
                os.mkdir(self._location)
        return self._location

    @location.setter
    def location(self, path_to_files: str):
        self._location = path_to_files

    @autocreate_directory_property
    def acquisitions_location(self) -> str:
        """str: location where acquisitions are stored"""
        return os.path.join(self.location, "acquisitions")

    @property
    def status(self) -> FileUploadStatus:
        """str: upload status based on the status of acquisitions"""
        child_status = set([f.status for f in self.acquisitions])
        if FileUploadStatus.UPLOADING in child_status:
            return FileUploadStatus.UPLOADING
        elif len(child_status) == 1 and FileUploadStatus.COMPLETE in child_status:
            return FileUploadStatus.COMPLETE
        else:
            return FileUploadStatus.WAITING

    def json(self):
        return {
            "id": self.id,
            "experiment_id": self.experiment_id,
            "name": self.name,
            "description": self.description,
            "location": self.location,
            "meta": self.meta,
            "created_at": self.created_at,
            "acquisitions": self.acquisitions,
        }

    def __repr__(self):
        return f"<Slide(id={self.id}, name={self.name})>"
