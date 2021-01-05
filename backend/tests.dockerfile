FROM python:3.8.7

LABEL maintainer="Anton Rau <anton.rau@gmail.com>"

ARG BACKEND_ENV=production

ENV PYTHONFAULTHANDLER=1 \
  PYTHONUNBUFFERED=1 \
  PYTHONHASHSEED=random \
  PIP_NO_CACHE_DIR=1 \
  PIP_DISABLE_PIP_VERSION_CHECK=on \
  PIP_DEFAULT_TIMEOUT=100 \
  POETRY_VERSION=1.1.4 \
  POETRY_VIRTUALENVS_CREATE=false \
  BACKEND_ENV=${BACKEND_ENV}
  PATH="${PATH}:/root/.poetry/bin" \
  PYTHONPATH=/app

# Install Poetry
RUN curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python

WORKDIR /app

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
COPY /histocat /app/pyproject.toml /app/
RUN poetry install --no-root

COPY histocat /app
RUN chmod +x /app/bin/tests-start.sh

# This will make the container wait, doing nothing, but alive
CMD ["bash", "-c", "while true; do sleep 1; done"]

# Afterwards you can exec a command /tests-start.sh in the live container, like:
# docker exec -it backend-tests /tests-start.sh
