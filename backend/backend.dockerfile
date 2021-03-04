FROM python:3.8.8

LABEL maintainer="Anton Rau <anton.rau@gmail.com>"

ARG BACKEND_ENV=production
ARG WORKERS_PER_CORE=0.5

ENV PYTHONFAULTHANDLER=1 \
  PYTHONUNBUFFERED=1 \
  PYTHONHASHSEED=random \
  PIP_NO_CACHE_DIR=1 \
  PIP_DISABLE_PIP_VERSION_CHECK=on \
  PIP_DEFAULT_TIMEOUT=100 \
  POETRY_VERSION=1.1.5 \
  POETRY_VIRTUALENVS_CREATE=false \
  BACKEND_ENV=${BACKEND_ENV} \
  PATH="${PATH}:/root/.poetry/bin" \
  PYTHONPATH=/app \
  WORKERS_PER_CORE=${WORKERS_PER_CORE}

# Install Poetry
RUN curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python

WORKDIR /app

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
COPY /poetry.lock /pyproject.toml /gunicorn_conf.py /app/
RUN poetry install --no-dev --no-root --extras backend

# Now copy in our code, and run it
COPY ./ /app
RUN chmod +x /app/bin/start.sh /app/bin/start-reload.sh

EXPOSE 80 5678

# Run the start script, it will check for an /app/prestart.sh script (e.g. for migrations)
# And then will start Gunicorn with Uvicorn
CMD ["/app/bin/start.sh"]

