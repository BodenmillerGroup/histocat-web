FROM python:3.7

WORKDIR /app

RUN pip install --no-cache pip-tools

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
COPY /app/requirements/prod.txt /app/requirements.txt
RUN pip-sync requirements.txt

# For development, Jupyter remote kernel, Hydrogen
# Using inside the container:
# jupyter notebook --ip=0.0.0.0 --allow-root
ARG env=prod
RUN bash -c "if [ $env == 'dev' ] ; then pip install jupyter ; fi"
EXPOSE 8888

ENV C_FORCE_ROOT=1

COPY ./app /app

ENV PYTHONPATH=/app

RUN chmod +x /app/worker-start.sh

CMD ["bash", "/app/worker-start.sh"]
