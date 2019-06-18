import logging
from datetime import timedelta, datetime

from sqlalchemy import Column, Integer, String, Interval, Boolean, BigInteger, LargeBinary, DateTime
from sqlalchemy.sql.functions import now

from app.db.base import Base

logger = logging.getLogger(__name__)


class Task(Base):
    """
    A *task* represents a computational job that can be submitted to
    computational resources for processing.
    Its state will be monitored while being processed and statistics will be
    collected.

    Warning
    -------
    This table is managed by GC3Pie. Don't modify it manually!
    """

    __tablename__ = "task"

    id: int = Column(Integer, primary_key=True, index=True)
    #: procssing state, e.g. ``"RUNNING"`` or ``"TERMINATED"``
    state: str = Column(String, index=True)
    #: name given by application
    name: str = Column(String, index=True)
    #: exit code
    exitcode: int = Column(Integer, index=True)
    #: total time of task
    time: timedelta = Column(Interval)
    #: total memory in MB of task
    memory: int = Column(Integer)
    #: total CPU time of task
    cpu_time: timedelta = Column(Interval)
    #: name of the corresponding Python object
    type: str = Column(String, index=True)
    #: whether the task is a collection of tasks
    is_collection: bool = Column(Boolean, index=True)
    #: ID of the parent task
    parent_id: int = Column(BigInteger, index=True)
    #: "Pickled" Python `gc3libs.Task` object
    data = Column(LargeBinary)
    #: ID of parent submission
    submission_id: int = Column(BigInteger, index=True)
    created_at: datetime = Column(DateTime, default=now(), nullable=False)
    updated_at: datetime = Column(DateTime, default=now(), onupdate=now(), nullable=False)

    def __repr__(self):
        return f'<Task(id={self.id}, name={self.name}, submission_id={self.submission_id})>'
