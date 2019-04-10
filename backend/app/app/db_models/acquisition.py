# TmLibrary - TissueMAPS library for distibuted image analysis routines.
# Copyright (C) 2016-2018 University of Zurich.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import logging
import os

from sqlalchemy import Column, String, Integer, Text, ForeignKey, UniqueConstraint
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship, backref, Session

from app.core.file_upload_status import FileUploadStatus
from app.core.utils import remove_location_upon_delete, autocreate_directory_property
from app.db_models import User
from app.db_models.base import DirectoryModel, CreatedAtMixin
from app.db_models.file import MicroscopeImageFile, MicroscopeMetadataFile

logger = logging.getLogger(__name__)

#: Format string for acquisition locations
ACQUISITION_LOCATION_FORMAT = 'acquisition_{id}'


@remove_location_upon_delete
class Acquisition(DirectoryModel, CreatedAtMixin):
    '''An *acquisition* contains all files belonging to one microscope image
    acquisition process. Note that in contrast to a *cycle*, an *acquisition*
    may contain more than one time point.

    The incentive to group files this way relates to the fact that most
    microscopes generate separate metadata files per *acquisition*.

    Attributes
    ----------
    microscope_image_files: List[tmlib.models.file.MicroscopeImageFile]
        microscope image files belonging to the acquisition
    microscope_metadata_files: List[tmlib.models.file.MicroscopeMetadataFile]
        microscope metadata files belonging to the acquisition
    channel_image_files: List[tmlib.models.file.ChannelImageFile]
        channel image files belonging to the acquisition
    '''

    __tablename__ = 'acquisition'
    __table_args__ = (UniqueConstraint('name', 'plate_id'),)

    #: str: name given by the user
    name = Column(String, index=True)

    #: str: description provided by the user
    description = Column(Text)

    #: int: ID of the parent plate
    plate_id = Column(
        Integer,
        ForeignKey('plate.id', onupdate='CASCADE', ondelete='CASCADE'),
        index=True
    )

    #: tmlib.models.plate.Plate: plate to which the acquisition belongs
    plate = relationship(
        'Plate',
        backref=backref('acquisition', cascade='all, delete-orphan')
    )

    def __init__(self, name: str, plate_id: int, description: str = ''):
        '''
        Parameters
        ----------
        name: str
            name of the acquisition
        plate_id: int
            ID of the parent :class:`Plate <tmlib.models.plate.Plate>`
        description: str, optional
            description of the acquisition
        '''
        # TODO: ensure that name is unique within plate
        self.name = name
        self.description = description
        self.plate_id = plate_id

    @hybrid_property
    def location(self) -> str:
        '''str: location were the acquisition content is stored'''
        if self._location is None:
            if self.id is None:
                raise AttributeError(
                    'Acquisition "%s" doesn\'t have an entry in the database yet. '
                    'Therefore, its location cannot be determined.' % self.name
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
    def microscope_images_location(self) -> str:
        '''str: location where microscope image files are stored'''
        return os.path.join(self.location, 'microscope_images')

    @autocreate_directory_property
    def microscope_metadata_location(self) -> str:
        '''str: location where microscope metadata files are stored'''
        return os.path.join(self.location, 'microscope_metadata')

    @property
    def status(self) -> FileUploadStatus:
        '''str: upload status based on the status of microscope files'''
        session = Session.object_session(self)
        img_files = session.query(MicroscopeImageFile.status). \
            filter_by(acquisition_id=self.id). \
            group_by(MicroscopeImageFile.status). \
            all()
        meta_files = session.query(MicroscopeMetadataFile.status). \
            filter_by(acquisition_id=self.id). \
            group_by(MicroscopeMetadataFile.status). \
            all()
        child_status = set([f.status for f in img_files]). \
            union([f.status for f in meta_files])
        if FileUploadStatus.UPLOADING in child_status:
            return FileUploadStatus.UPLOADING
        elif len(child_status) == 1 and FileUploadStatus.COMPLETE in child_status:
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

    def belongs_to(self, user: User) -> bool:
        '''Determines whether the acquisition belongs to a given `user`.

        Parameters
        ----------
        user: tmlib.user.User
            `TissueMAPS` user

        Returns
        -------
        bool
            whether acquisition belongs to `user`
        '''
        return self.plate.belongs_to(user)

    def __repr__(self):
        return f'<Acquisition(id={self.id}, name={self.name})>'
