FROM python:3.7

ENV PIP_DISABLE_PIP_VERSION_CHECK=on

# Install Poetry
RUN pip install --no-cache poetry==0.12.12 && poetry config settings.virtualenvs.create false

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
WORKDIR /app
COPY /app/poetry.lock /app/pyproject.toml /app/
RUN poetry install --no-interaction --no-dev

ENV C_FORCE_ROOT=1

COPY ./app /app

ENV PYTHONPATH=/app

RUN chmod +x /app/worker-start.sh

CMD ["bash", "/app/worker-start.sh"]
