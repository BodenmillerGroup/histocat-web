import logging
import os

import dramatiq
import emails
from dramatiq.brokers.rabbitmq import RabbitmqBroker
from emails.template import JinjaTemplate

from app.core import config
from app.db.session import db_session
from app.io.mcd import import_mcd
from app.io.ome_tiff_loader import OmeTiffLoader
from app.io.text_loader import TextLoader

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

        import pydevd_pycharm
        # TODO: Don't forget to modify IP address!!
        pydevd_pycharm.settrace('130.60.106.36', port=5679, stdoutToServer=True, stderrToServer=True)
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


@dramatiq.actor(queue_name='processing')
def import_slide(uri: str, experiment_id: int):
    logger.info(f'Importing slide into experiment [{experiment_id}] from {uri}')

    filename, file_extension = os.path.splitext(uri)
    file_extension = file_extension.lower()

    if file_extension == ".mcd":
        import_mcd(db_session, uri, experiment_id)
    elif file_extension == ".tiff" or file_extension == ".tif":
        if filename.endswith(".ome"):
            OmeTiffLoader.load(db_session, uri, experiment_id)
    elif file_extension == ".txt":
        TextLoader.load(db_session, uri, experiment_id)

    os.remove(uri)


@dramatiq.actor(queue_name='processing')
def prepare_dataset(dataset_id: int):
    logger.info(f'Preparing dataset [{dataset_id}]...')
    dataset = dataset_crud.get(db_session, id=dataset_id)
    logger.info(dataset.meta["input"])
