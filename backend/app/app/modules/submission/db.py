import logging
from datetime import datetime

from sqlalchemy import Column, Integer, String, BigInteger, ForeignKey, DateTime
from sqlalchemy.orm import relationship, backref
from sqlalchemy.sql.functions import now

from app.db.base import Base

logger = logging.getLogger(__name__)


class Submission(Base):
    """
    A *submission* handles the processing of a computational *task* on a cluster.
    """

    __tablename__ = 'submission'

    id: int = Column(BigInteger, primary_key=True, autoincrement=True, index=True)
    #: name of the program that submitted the tasks
    program: str = Column(String, index=True)
    #: ID of the parent experiment
    experiment_id: int = Column(
        Integer,
        ForeignKey("experiment.id", ondelete="CASCADE"),
        index=True,
    )
    #: ID of the submitting user
    user_id: int = Column(
        Integer,
        ForeignKey("user.id", ondelete="CASCADE"),
        index=True,
    )
    # TODO: make top_task_id a foreign key and create a relationship
    #: ID of the top task in the submitted collection of tasks
    top_task_id: int = Column(BigInteger, index=True)
    created_at: datetime = Column(DateTime, default=now(), nullable=False)

    experiment = relationship(
        'Experiment',
        backref=backref('submissions', cascade='all, delete-orphan')
    )

    user = relationship(
        'User',
        backref=backref('submissions', cascade='all, delete-orphan')
    )

    def __init__(self, experiment_id: int, program: str, user_id: int):
        self.experiment_id = experiment_id
        self.program = program
        self.user_id = user_id

    def __repr__(self):
        return f'<Submission(id={self.id}, experiment_id={self.experiment_id}, program={self.program}, user_id={self.user_id})>'
