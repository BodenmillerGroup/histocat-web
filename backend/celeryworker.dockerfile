FROM python:3.7.3

LABEL maintainer="anton.rau@gmail.com"

ENV PYTHONFAULTHANDLER=1 \
  PYTHONUNBUFFERED=1 \
  PYTHONHASHSEED=random \
  PIP_NO_CACHE_DIR=off \
  PIP_DISABLE_PIP_VERSION_CHECK=on \
  PIP_DEFAULT_TIMEOUT=100 \
  POETRY_VERSION=0.12.16

# Install Poetry
RUN pip install --no-cache "poetry==$POETRY_VERSION" && poetry config settings.virtualenvs.create false

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
