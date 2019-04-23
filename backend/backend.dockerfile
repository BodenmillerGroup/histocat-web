FROM tiangolo/uvicorn-gunicorn:python3.7

WORKDIR /app

RUN pip install --no-cache pip-tools

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
COPY /app/requirements/dev.txt /app/requirements.txt
RUN pip-sync requirements.txt

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
