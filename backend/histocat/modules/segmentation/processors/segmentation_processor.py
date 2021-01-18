import tensorflow as tf
from fastapi import HTTPException
from sqlalchemy.orm import Session
from starlette import status

from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.core.utils import timeit
from histocat.modules.model import service as model_service
from histocat.modules.project import service as project_service
from histocat.modules.segmentation.dto import SegmentationSubmissionDto
from histocat.modules.segmentation.processors import acquisition_processor


@timeit
def process_segmentation(db: Session, project_id: int, params: SegmentationSubmissionDto):
    """
    Segmentation processing
    """

    project = project_service.get(db, id=project_id)
    if not project:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Project id:{project_id} not found")

    model = model_service.get(db, id=params.model_id)
    if not model:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Model id:{params.model_id} not found")

    # disable GPU
    tf.config.set_visible_devices([], "GPU")

    # load UNET model
    keras_model = tf.keras.models.load_model(model.location)

    for acquisition_id in params.acquisition_ids:
        acquisition_processor.process_acquisition(db, acquisition_id=acquisition_id, params=params, model=keras_model)

    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(project_id, "segmentation_ready", None))
