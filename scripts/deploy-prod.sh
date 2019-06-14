#! /usr/bin/env bash

# Exit in case of error
set -e

BACKEND_ENV=production \
FRONTEND_ENV=production \
VUE_APP_ENV=production \
docker-compose \
-f docker/shared.yml \
-f docker/prod.yml \
config > docker-stack.yml

docker-compose -f docker-stack.yml up -d
