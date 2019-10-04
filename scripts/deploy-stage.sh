#! /usr/bin/env bash

# Exit in case of error
set -e

BACKEND_ENV=stage \
FRONTEND_ENV=stage \
VUE_APP_ENV=stage \
docker-compose \
-f docker/shared.yml \
-f docker/stage.yml \
config > docker-stack.yml

docker-compose -f docker-stack.yml up -d
