from sqlalchemy import Column, Integer, String, Text, DateTime
from sqlalchemy.sql.functions import now

from app.db_models.base import Base

#: Format string for experiment locations.
EXPERIMENT_LOCATION_FORMAT = 'experiment_{id}'


class Experiment(Base):
    __tablename__ = 'experiment'

    id = Column(Integer, primary_key=True, autoincrement=True, index=True)
    name = Column(String, index=True, nullable=False)
    description = Column(Text, nullable=True)
    #: absolute path to the directory where experiments are located
    root_directory = Column(String)
    created_at = Column(DateTime, default=now())

    def __repr__(self):
        return f'<Experiment(id={self.id}, name={self.name}, description={self.description},' \
            f' root_directory={self.root_directory})>'
