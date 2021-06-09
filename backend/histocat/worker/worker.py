import logging
import os
import shutil

import dramatiq
import emails
from dramatiq.brokers.rabbitmq import RabbitmqBroker
from emails.template import JinjaTemplate
from imctools.io.utils import MCD_FILENDING, ZIP_FILENDING

import histocat.core.init_db  # noqa
from histocat.config import config
from histocat.core.errors import DataImportError
from histocat.core.pipeline.dto import PipelineProcessDto
from histocat.core.segmentation.dto import SegmentationSubmissionDto
from histocat.core.session import db_session
from histocat.worker.io import mcd, model, zip
from histocat.worker.pipeline.processors import pipeline_processor
from histocat.worker.segmentation.processors import segmentation_processor

rabbitmq_broker = RabbitmqBroker(host="rabbitmq", connection_attempts=10)
dramatiq.set_broker(rabbitmq_broker)

logger = logging.getLogger(__name__)

if os.environ.get("BACKEND_ENV") == "development":
    try:
        # VS Code Debugging
        # Allow other computers to attach to ptvsd at this IP address and port.
        import ptvsd

        ptvsd.enable_attach(address=("0.0.0.0", 5688), redirect_output=True)
        pass

        # PyCharm Debugging
        # TODO: Don't forget to modify IP address!!
        # import pydevd_pycharm
        # pydevd_pycharm.settrace('130.60.106.48', port=5679, stdoutToServer=True, stderrToServer=True, suspend=False)

    except Exception as e:
        logger.error(e)


@dramatiq.actor(queue_name="default")
def test_worker(word: str):
    logger.info(f"Testing worker: [{word}]")


@dramatiq.actor(queue_name="email")
def send_email(email_to: str, subject_template="", html_template="", environment={}):
    logger.info(f"Sending email to: [{email_to}]")
    assert config.EMAILS_ENABLED, "no provided configuration for email variables"
    message = emails.Message(
        subject=JinjaTemplate(subject_template),
        html=JinjaTemplate(html_template),
        mail_from=(config.EMAILS_FROM_NAME, config.EMAILS_FROM),
    )
    smtp_options = {"host": config.SMTP_HOST, "port": config.SMTP_PORT}
    if config.SMTP_TLS:
        smtp_options["tls"] = True
    if config.SMTP_USER:
        smtp_options["user"] = config.SMTP_USER
    if config.SMTP_PASSWORD:
        smtp_options["password"] = config.SMTP_PASSWORD
    response = message.send(to=email_to, render=environment, smtp=smtp_options)
    logging.info(f"Send email result: {response}")


@dramatiq.actor(queue_name="import", max_retries=0, time_limit=1000 * 60 * 60 * 10)  # 10 hours time limit
def import_slide(uri: str, project_id: int):
    logger.info(f"Importing slide into project [{project_id}] from {uri}")

    path = os.path.dirname(os.path.abspath(uri))
    filename, file_extension = os.path.splitext(uri)
    file_extension = file_extension.lower()

    try:
        if file_extension == MCD_FILENDING:
            mcd.import_mcd(db_session, uri, project_id)
        elif file_extension == ZIP_FILENDING:
            zip.import_slide(db_session, uri, project_id)
    except DataImportError as error:
        logger.warning(error)
    finally:
        shutil.rmtree(path)


@dramatiq.actor(queue_name="import", max_retries=0, time_limit=1000 * 60 * 60 * 10)  # 10 hours time limit
def import_dataset(type: str, masks_folder: str, regionprops_folder: str, intensities_folder: str, uri: str, project_id: int):
    logger.info(f"Importing dataset into project [{project_id}] from {uri}")

    path = os.path.dirname(os.path.abspath(uri))
    filename, file_extension = os.path.splitext(uri)
    file_extension = file_extension.lower()

    try:
        if file_extension == ZIP_FILENDING:
            zip.import_dataset(db_session, type, masks_folder, regionprops_folder, intensities_folder, uri, project_id)
    except DataImportError as error:
        logger.warning(error)
    finally:
        shutil.rmtree(path)


@dramatiq.actor(queue_name="process", max_retries=0, time_limit=1000 * 60 * 60 * 10)  # 10 hours time limit
def process_pipeline(payload: str):
    params: PipelineProcessDto = PipelineProcessDto.parse_raw(payload)
    logger.info(f"Processing pipeline for acquisitions {params.acquisition_ids} from dataset [{params.dataset_id}]")
    try:
        pipeline_processor.process_pipeline(db_session, params)
    except Exception as error:
        logger.warning(error)
    finally:
        pass


@dramatiq.actor(queue_name="import", max_retries=0, time_limit=1000 * 60 * 60 * 10)  # 10 hours time limit
def import_model(uri: str, model_id: int):
    logger.info(f"Importing model [{model_id}] from {uri}")

    path = os.path.dirname(os.path.abspath(uri))
    filename, file_extension = os.path.splitext(uri)

    try:
        if file_extension.lower() == ZIP_FILENDING:
            model.import_model_zip(db_session, uri, model_id)
    except DataImportError as error:
        logger.warning(error)
    finally:
        shutil.rmtree(path)


@dramatiq.actor(queue_name="process", max_retries=0, time_limit=1000 * 60 * 60 * 10)  # 10 hours time limit
def process_segmentation(project_id: int, payload: str):
    params: SegmentationSubmissionDto = SegmentationSubmissionDto.parse_raw(payload)
    logger.info(f"Processing segmentation for acquisitions {params.acquisition_ids} with model [{params.model_id}]")
    try:
        segmentation_processor.process_segmentation(db_session, project_id=project_id, params=params)
    except Exception as error:
        logger.warning(error)
    finally:
        pass
