FROM python:3.7.6

LABEL maintainer="Anton Rau <anton.rau@gmail.com>"

ARG BACKEND_ENV=production
ARG WORKERS_PER_CORE=0.5

ENV PYTHONFAULTHANDLER=1 \
  PYTHONUNBUFFERED=1 \
  PYTHONHASHSEED=random \
  PIP_NO_CACHE_DIR=off \
  PIP_DISABLE_PIP_VERSION_CHECK=on \
  PIP_DEFAULT_TIMEOUT=100 \
  POETRY_VERSION=1.0.2 \
  BACKEND_ENV=${BACKEND_ENV} \
  WORKERS_PER_CORE=${WORKERS_PER_CORE}

# Install Poetry
RUN pip install --no-cache "poetry==$POETRY_VERSION" && poetry config virtualenvs.create false

WORKDIR /app

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
COPY /app/poetry.lock /app/pyproject.toml /app/gunicorn_conf.py /app/
RUN poetry install --no-dev --extras backend

# Now copy in our code, and run it
COPY ./app /app
RUN chmod +x /app/bin/start.sh /app/bin/start-reload.sh

ENV PYTHONPATH=/app

EXPOSE 80 5678

# Run the start script, it will check for an /app/prestart.sh script (e.g. for migrations)
# And then will start Gunicorn with Uvicorn
CMD ["/app/bin/start.sh"]

