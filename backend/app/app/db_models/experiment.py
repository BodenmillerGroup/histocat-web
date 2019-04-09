from datetime import datetime

from sqlalchemy import Boolean, Column, Integer, String, Text, ForeignKey, DateTime
from sqlalchemy.orm import relationship

from app.db.base_class import Base


class Experiment(Base):
    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, index=True)
    description = Column(Text)
    #: absolute path to the directory where experiments are located
    root_directory = Column(String)
    #: int: ID of the owner
    owner_id = Column(Integer, ForeignKey('user.id'), index=True)
    #: db_models.user.User: user that owns the experiment
    owner = relationship('User', back_populates='experiment')
    created_at = Column(DateTime, default=datetime.utcnow)

    def __repr__(self):
        return f'<Experiment(id={self.id}, name={self.name}, description={self.description},' \
            f' root_directory={self.root_directory}, owner_id={self.owner_id})>'
