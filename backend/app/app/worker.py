import logging
import os
import shutil

import dramatiq
import emails
from dramatiq.brokers.rabbitmq import RabbitmqBroker
from emails.template import JinjaTemplate
from imctools.scripts.convertfolder2imcfolder import MCD_FILENDING, ZIP_FILENDING
from typing import List

from app.core import config
from app.core.errors import SlideImportError
from app.db.session import db_session
from app.io import mcd
from app.io import zip
from app.modules.analysis.processors import tsne, umap, phenograph

import app.db.init_db  # noqa

rabbitmq_broker = RabbitmqBroker(host="rabbitmq", connection_attempts=10)
dramatiq.set_broker(rabbitmq_broker)

logger = logging.getLogger(__name__)

if os.environ.get("BACKEND_ENV") == "development":
    try:
        # VS Code Debugging
        # Allow other computers to attach to ptvsd at this IP address and port.
        # import ptvsd
        # ptvsd.enable_attach(address=('0.0.0.0', 5688), redirect_output=True)

        # PyCharm Debugging
        # TODO: Don't forget to modify IP address!!
        # import pydevd_pycharm
        # pydevd_pycharm.settrace('130.60.106.31', port=5679, stdoutToServer=True, stderrToServer=True, suspend=False)

        pass
    except Exception as e:
        logger.error(e)


@dramatiq.actor(queue_name='default')
def test_worker(word: str):
    logger.info(f'Testing worker: [{word}]')


@dramatiq.actor(queue_name='email')
def send_email(email_to: str, subject_template="", html_template="", environment={}):
    logger.info(f'Sending email to: [{email_to}]')
    assert config.EMAILS_ENABLED, "no provided configuration for email variables"
    message = emails.Message(
        subject=JinjaTemplate(subject_template),
        html=JinjaTemplate(html_template),
        mail_from=(config.EMAILS_FROM_NAME, config.EMAILS_FROM_EMAIL),
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


@dramatiq.actor(queue_name='import', max_retries=0)
def import_data(uri: str, experiment_id: int, user_id: int):
    logger.info(f'Importing data into experiment [{experiment_id}] from {uri}')

    path = os.path.dirname(os.path.abspath(uri))
    filename, file_extension = os.path.splitext(uri)
    file_extension = file_extension.lower()

    try:
        if file_extension == MCD_FILENDING:
            mcd.import_mcd(db_session, uri, experiment_id, user_id)
        elif file_extension == ZIP_FILENDING:
            zip.import_zip(db_session, uri, experiment_id, user_id)
    except SlideImportError as error:
        logger.warning(error)
    finally:
        shutil.rmtree(path)


@dramatiq.actor(queue_name='process', max_retries=0, time_limit=1000 * 60 * 60 * 10)
def process_tsne(
    dataset_id: int,
    acquisition_ids: List[int],
    n_components: int,
    perplexity: int,
    learning_rate: int,
    iterations: int,
    theta: float,
    init: str,
    markers: List[str],
):
    logger.info(f'Processing t-SNE for acquisitions {acquisition_ids} from dataset [{dataset_id}]')
    try:
        tsne.process_tsne(db_session, dataset_id, acquisition_ids, n_components, perplexity, learning_rate, iterations, theta, init, markers)
    except Exception as error:
        logger.warning(error)
    finally:
        pass


@dramatiq.actor(queue_name='process', max_retries=0, time_limit=1000 * 60 * 60 * 10)
def process_umap(
    dataset_id: int,
    acquisition_ids: List[int],
    n_components: int,
    n_neighbors: int,
    metric: str,
    min_dist: float,
    markers: List[str],
):
    logger.info(f'Processing UMAP for acquisitions {acquisition_ids} from dataset [{dataset_id}]')
    try:
        umap.process_umap(db_session, dataset_id, acquisition_ids, n_components, n_neighbors, metric, min_dist, markers)
    except Exception as error:
        logger.warning(error)
    finally:
        pass


@dramatiq.actor(queue_name='process', max_retries=0, time_limit=1000 * 60 * 60 * 10)
def process_phenograph(
    dataset_id: int,
    acquisition_ids: List[int],
    markers: List[str],
):
    logger.info(f'Processing PhenoGraph for acquisitions {acquisition_ids} from dataset [{dataset_id}]')
    try:
        phenograph.process_phenograph(db_session, dataset_id, acquisition_ids, markers)
    except Exception as error:
        logger.warning(error)
    finally:
        pass
