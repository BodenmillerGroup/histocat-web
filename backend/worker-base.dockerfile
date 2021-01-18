FROM ubuntu:20.04

LABEL maintainer="Anton Rau <anton.rau@gmail.com>"

ENV PYTHONFAULTHANDLER=1 \
  PYTHONUNBUFFERED=1 \
  PYTHONHASHSEED=random \
  PIP_NO_CACHE_DIR=1 \
  PIP_DISABLE_PIP_VERSION_CHECK=on \
  PIP_DEFAULT_TIMEOUT=100 \
  POETRY_VERSION=1.1.4 \
  POETRY_VIRTUALENVS_CREATE=false \
  PATH="${PATH}:/root/.poetry/bin" \
  PYTHONPATH=/app

# update system
RUN apt-get update && apt-get upgrade -y

# install packages
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt install -y \
    software-properties-common \
    apt-utils \
    build-essential \
    default-libmysqlclient-dev \
    git \
    libgtk-3-dev \
    libnotify-dev \
    libsdl2-dev \
    locales \
    curl \
    openjdk-11-jdk-headless \
    python3.8-dev \
    python3.8-distutils \
    python3-pip \
    wget && apt clean

# set java paths
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
ENV PATH=$PATH:/home/ubuntu/.local/bin

# Install Poetry
RUN curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python3

WORKDIR /app

# Set the locale
RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
COPY /poetry.lock /pyproject.toml /app/
RUN poetry install --no-dev --no-root --extras worker
