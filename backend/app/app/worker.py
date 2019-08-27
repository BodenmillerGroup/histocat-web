import datetime
import logging
import os
import shutil

import dramatiq
import emails
from dramatiq.brokers.rabbitmq import RabbitmqBroker
from emails.template import JinjaTemplate
from imctools.scripts.convertfolder2imcfolder import MCD_FILENDING, ZIP_FILENDING

from app.core import config
from app.core.errors import SlideImportError
from app.core.notifier import Message
from app.core.redis_manager import redis_manager, UPDATES_CHANNEL_NAME
from app.db.session import db_session
from app.io import dataset, mcd
from app.io import zip

from app.modules.dataset import crud as dataset_crud

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
        # import pydevd_pycharm
        # TODO: Don't forget to modify IP address!!
        # pydevd_pycharm.settrace('130.60.106.25', port=5679, stdoutToServer=True, stderrToServer=True, suspend=False)

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
def import_slide(uri: str, experiment_id: int):
    logger.info(f'Importing slide into experiment [{experiment_id}] from {uri}')

    path = os.path.dirname(os.path.abspath(uri))
    filename, file_extension = os.path.splitext(uri)
    file_extension = file_extension.lower()

    try:
        if file_extension == MCD_FILENDING:
            mcd.import_mcd(db_session, uri, experiment_id)
        elif file_extension == ZIP_FILENDING:
            zip.import_zip(db_session, uri, experiment_id)
    except SlideImportError as error:
        logger.warn(error)
    finally:
        shutil.rmtree(path)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(experiment_id, "slide_imported"))


@dramatiq.actor(queue_name='import', max_retries=0)
def prepare_dataset(dataset_id: int):
    logger.info(f'Preparing dataset [{dataset_id}].')
    item = dataset_crud.get(db_session, id=dataset_id)

    if not item:
        logger.error(f'Dataset [{dataset_id}] not found')
        return

    try:
        dataset.prepare_dataset(db_session, item)
        item.status = 'imported'
        db_session.add(item)
        db_session.commit()
    except Exception as error:
        if not item.errors:
            item.errors = dict()
        item.errors[str(datetime.datetime.now())] = str(error)
        db_session.add(item)
        db_session.commit()
        logger.error(error)


@dramatiq.actor(queue_name='import', max_retries=0)
def import_dataset(uri: str, user_id: int, experiment_id: int):
    logger.info(f'Importing dataset into experiment [{experiment_id}] from {uri} by user {user_id}')

    path = os.path.dirname(os.path.abspath(uri))
    filename, file_extension = os.path.splitext(uri)
    file_extension = file_extension.lower()

    try:
        if file_extension == ZIP_FILENDING:
            dataset.import_zip(db_session, uri, user_id, experiment_id)
    except SlideImportError as error:
        logger.warn(error)
    finally:
        shutil.rmtree(path)
