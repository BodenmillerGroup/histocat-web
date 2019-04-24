FROM tiangolo/uvicorn-gunicorn:python3.7

ENV PIP_DISABLE_PIP_VERSION_CHECK=on

# Install Poetry
RUN pip install --no-cache poetry==0.12.12 && poetry config settings.virtualenvs.create false

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
WORKDIR /app
COPY /app/poetry.lock /app/pyproject.toml /app/
RUN poetry install --no-interaction --no-dev

# For development, Jupyter remote kernel, Hydrogen
# Using inside the container:
# jupyter notebook --ip=0.0.0.0 --allow-root
ARG env=prod
RUN bash -c "if [ $env == 'dev' ] ; then pip install jupyter ; fi"
EXPOSE 8888

# Now copy in our code, and run it
COPY ./app /app

ENV PYTHONPATH=/app

EXPOSE 80 5678
