#! /usr/bin/env bash

# Exit in case of error
set -e

if [ $(uname -s) = "Linux" ]; then
    echo "Remove __pycache__ files"
    sudo find . -type d -name __pycache__ -exec rm -r {} \+
fi

JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64 \
docker build \
-f ./backend/worker-base.dockerfile \
-t plankter/histocat-worker-base \
./backend
docker push plankter/histocat-worker-base
