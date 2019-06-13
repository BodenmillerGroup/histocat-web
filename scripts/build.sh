#! /usr/bin/env bash

# Exit in case of error
set -e

if [ $(uname -s) = "Linux" ]; then
    echo "Remove __pycache__ files"
    sudo find . -type d -name __pycache__ -exec rm -r {} \+
fi

TAG=${TAG-latest} \
BACKEND_ENV=${BACKEND_ENV-production} \
FRONTEND_ENV=${FRONTEND_ENV-production} \
docker-compose \
-f docker/deploy.build.yml \
-f docker/deploy.images.yml \
config > docker-stack.yml

docker-compose -f docker-stack.yml build
