FROM python:3.7.6

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
RUN poetry install --no-dev

COPY ./app /app
RUN chmod +x /app/bin/worker-start.sh

ENV PYTHONPATH=/app

EXPOSE 5688

CMD ["/app/bin/worker-start.sh"]
