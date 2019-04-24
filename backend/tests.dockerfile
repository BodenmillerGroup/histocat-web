FROM python:3.7

ENV PIP_DISABLE_PIP_VERSION_CHECK=on

# Install Poetry
RUN pip install --no-cache poetry==0.12.12 && poetry config settings.virtualenvs.create false

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
