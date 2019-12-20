FROM python:3.7.5

LABEL maintainer="Anton Rau <anton.rau@gmail.com>"

ARG BACKEND_ENV=production

ENV PYTHONFAULTHANDLER=1 \
  PYTHONUNBUFFERED=1 \
  PYTHONHASHSEED=random \
  PIP_NO_CACHE_DIR=off \
  PIP_DISABLE_PIP_VERSION_CHECK=on \
  PIP_DEFAULT_TIMEOUT=100 \
  POETRY_VERSION=1.0.0 \
  BACKEND_ENV=${BACKEND_ENV}

# Install Poetry
RUN pip install --no-cache "poetry==$POETRY_VERSION" && poetry config virtualenvs.create false

WORKDIR /app

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
COPY /app/poetry.lock /app/pyproject.toml /app/
RUN poetry install

COPY ./app /app
RUN chmod +x /app/bin/tests-start.sh

ENV PYTHONPATH=/app

# This will make the container wait, doing nothing, but alive
CMD ["bash", "-c", "while true; do sleep 1; done"]

# Afterwards you can exec a command /tests-start.sh in the live container, like:
# docker exec -it backend-tests /tests-start.sh
