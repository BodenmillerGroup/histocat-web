#! /usr/bin/env bash

# Exit in case of error
set -e

DOMAIN=localhost \
BACKEND_ENV=development \
FRONTEND_ENV=development \
VUE_APP_ENV=development \
docker-compose \
-f docker/shared.yml \
-f docker/dev.yml \
config > docker-stack.yml

docker-compose -f docker-stack.yml up -d
