FROM python:3.7.3

LABEL maintainer="anton.rau@gmail.com"

ARG BACKEND_ENV=production

ENV PYTHONFAULTHANDLER=1 \
  PYTHONUNBUFFERED=1 \
  PYTHONHASHSEED=random \
  PIP_NO_CACHE_DIR=off \
  PIP_DISABLE_PIP_VERSION_CHECK=on \
  PIP_DEFAULT_TIMEOUT=100 \
  POETRY_VERSION=0.12.16 \
  BACKEND_ENV=${BACKEND_ENV}

# Install Poetry
RUN pip install --no-cache "poetry==$POETRY_VERSION" && poetry config settings.virtualenvs.create false

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
WORKDIR /app
COPY /app/poetry.lock /app/pyproject.toml /app/
RUN poetry install --no-interaction

COPY ./app /app

ENV PYTHONPATH=/app

RUN chmod +x /app/tests-start.sh

# This will make the container wait, doing nothing, but alive
CMD ["bash", "-c", "while true; do sleep 1; done"]

# Afterwards you can exec a command /tests-start.sh in the live container, like:
# docker exec -it backend-tests /tests-start.sh
