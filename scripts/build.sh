#! /usr/bin/env bash

# Exit in case of error
set -e

TAG=${TAG} \
FRONTEND_ENV=${FRONTEND_ENV-production} \
docker-compose \
-f docker/deploy.build.yml \
-f docker/deploy.images.yml \
config > docker-stack.yml

docker-compose -f docker-stack.yml build
