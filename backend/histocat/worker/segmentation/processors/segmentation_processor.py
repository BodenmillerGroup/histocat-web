import tensorflow as tf
from sqlalchemy.orm import Session

from histocat.core.dataset import service as dataset_service
from histocat.core.dataset.dto import DatasetCreateDto, DatasetUpdateDto
from histocat.core.errors import SegmentationError
from histocat.core.model import service as model_service
from histocat.core.notifier import Message
from histocat.core.project import service as project_service
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.core.segmentation.dto import SegmentationSubmissionDto
from histocat.core.utils import timeit
from histocat.worker.segmentation.processors import (
    acquisition_processor,
    dataset_importer,
)


@timeit
def process_segmentation(db: Session, project_id: int, params: SegmentationSubmissionDto):
    """
    Segmentation processing
    """

    project = project_service.get(db, id=project_id)
    if not project:
        raise SegmentationError(f"Project id:{project_id} not found")

    model = model_service.get(db, id=params.model_id)
    if not model:
        raise SegmentationError(f"Model id:{params.model_id} not found")

    # disable GPU
    tf.config.set_visible_devices([], "GPU")

    # Set CPU cores limit
    # tf.config.threading.set_intra_op_parallelism_threads(4)
    # tf.config.threading.set_inter_op_parallelism_threads(4)

    # load UNET model
    keras_model = tf.keras.models.load_model(model.location, compile=False)

    dataset_params = DatasetCreateDto(project_id=project_id, status="pending")
    if params.dataset_name:
        dataset_params.name = params.dataset_name
    if params.dataset_description:
        dataset_params.description = params.dataset_description
    dataset = dataset_service.create(db, params=dataset_params)

    segmentation_data = []
    for acquisition_id in params.acquisition_ids:
        ac_output = acquisition_processor.process_acquisition(
            db, acquisition_id=acquisition_id, params=params, model=keras_model, dataset=dataset
        )
        segmentation_data.append(ac_output)

    meta = dataset_importer.import_dataset(db, dataset=dataset, segmentation_data=segmentation_data)
    meta["params"] = params.dict()

    dataset_service.update(db, item=dataset, params=DatasetUpdateDto(status="ready", meta=meta))

    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(project_id, "segmentation_ready", None))
