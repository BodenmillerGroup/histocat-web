FROM tiangolo/uvicorn-gunicorn:python3.7

LABEL maintainer="anton.rau@gmail.com"

ENV PYTHONFAULTHANDLER=1 \
  PYTHONUNBUFFERED=1 \
  PYTHONHASHSEED=random \
  PIP_NO_CACHE_DIR=off \
  PIP_DISABLE_PIP_VERSION_CHECK=on \
  PIP_DEFAULT_TIMEOUT=100 \
  POETRY_VERSION=0.12.16 \
  ENV=dev

# Install Poetry
RUN pip install --no-cache "poetry==$POETRY_VERSION" && poetry config settings.virtualenvs.create false

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
WORKDIR /app
COPY /app/poetry.lock /app/pyproject.toml /app/
RUN poetry install --no-interaction

# Now copy in our code, and run it
COPY ./app /app

ENV PYTHONPATH=/app

EXPOSE 80 5678
