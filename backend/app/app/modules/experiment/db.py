import logging
import os

from sqlalchemy import Column, String, Text
from sqlalchemy.orm import relationship

from app.core.utils import autocreate_directory_property, remove_location_upon_delete
from app.db.base import CreatedAtMixin, DirectoryModel, MetaMixin

logger = logging.getLogger(__name__)

#: Format string for experiment locations.
EXPERIMENT_LOCATION_FORMAT = "experiment_{id}"


@remove_location_upon_delete
class Experiment(DirectoryModel, MetaMixin, CreatedAtMixin):
    """An *experiment* is the main organizational unit of `HistoCAT`.
    It represents a set of images and associated data.

    Attributes
    ----------
    slides: List[app.db_models.slide.Slide]
        slides belonging to the experiment
    """

    __tablename__ = "experiment"

    #: str: name given by the user
    name = Column(String, index=True)

    #: str: description provided by the user
    description = Column(Text)

    #: absolute path to the directory where experiments are located
    root_directory = Column(String)

    slides = relationship("Slide", back_populates="experiment")

    def __init__(
        self, name: str, root_directory: str, description: str = "", meta: dict = None
    ):
        """
        Parameters
        ----------
        name: str
            name of the experiment
        root_directory: str
            absolute path to root directory on disk where experiment should
            be created in
        description: str, optional
            description of the experimental setup
        meta: dict, optional
            meta data of the experiment
        """
        self.name = name
        self.description = description
        self.root_directory = root_directory
        self.meta = meta

    @autocreate_directory_property
    def location(self):
        """str: location of the experiment,
        e.g. absolute path to a directory on disk
        """
        if self.id is None:
            raise AttributeError(
                f'Experiment "{self.name}" doesn\'t have an entry in the database yet. Therefore, its location cannot be determined.'
            )
        return os.path.join(
            os.path.expandvars(self.root_directory),
            EXPERIMENT_LOCATION_FORMAT.format(id=self.id),
        )

    @autocreate_directory_property
    def slides_location(self):
        """str: location where slides data are stored"""
        return os.path.join(self.location, "slides")

    def json(self):
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "root_directory": self.root_directory,
            "location": self.location,
            "meta": self.meta,
            "created_at": self.created_at,
            "slides": self.slides,
        }

    def __repr__(self):
        return "<Experiment(id=%r, name=%r)>" % (self.id, self.name)
