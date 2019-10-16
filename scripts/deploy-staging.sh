#! /usr/bin/env bash

# Exit in case of error
set -e

BACKEND_ENV=staging \
FRONTEND_ENV=staging \
VUE_APP_ENV=staging \
docker-compose \
-f .deploy/shared.yml \
-f .deploy/staging.yml \
config > docker-stack.yml

docker-compose -f docker-stack.yml up -d
