import logging
import os

import dramatiq
import emails
from dramatiq.brokers.rabbitmq import RabbitmqBroker
from emails.template import JinjaTemplate

import app.db.init_db  # noqa
from app.core import config
from app.db.session import db_session
from app.io.mcd_loader import McdLoader
from app.io.ome_tiff_loader import OmeTiffLoader

rabbitmq_broker = RabbitmqBroker(host="rabbitmq", connection_attempts=10)
dramatiq.set_broker(rabbitmq_broker)

logger = logging.getLogger(__name__)


@dramatiq.actor(queue_name='default')
def test_worker(word: str):
    logger.info(f'Testing worker: [{word}]')


@dramatiq.actor(queue_name='default')
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


@dramatiq.actor(queue_name='default')
def import_slide(uri: str, experiment_id: int):
    logger.info(f'Importing slide into experiment [{experiment_id}] from {uri}')

    filename, file_extension = os.path.splitext(uri)
    file_extension = file_extension.lower()

    if file_extension == ".mcd":
        McdLoader.load(db_session, uri, experiment_id)
    elif file_extension == ".tiff" or file_extension == ".tif":
        if filename.endswith(".ome"):
            OmeTiffLoader.load(db_session, uri, experiment_id)
    elif file_extension == ".txt":
        pass

    os.remove(uri)
