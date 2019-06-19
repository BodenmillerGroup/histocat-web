import logging

import dramatiq
import emails
from dramatiq.brokers.rabbitmq import RabbitmqBroker
from emails.template import JinjaTemplate

from app.core import config

rabbitmq_broker = RabbitmqBroker(host="rabbitmq", connection_attempts=10)
dramatiq.set_broker(rabbitmq_broker)

logger = logging.getLogger(__name__)


@dramatiq.actor(queue_name='default')
def test_worker(word: str):
    logger.info(f'Testing worker: [{word}]')


@dramatiq.actor(queue_name='default')
def import_slide(experiment_id: int, uri: str):
    logger.info(f'Importing slide: [{experiment_id}, {uri}]')


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
