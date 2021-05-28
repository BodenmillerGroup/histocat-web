#! /usr/bin/env bash

# Exit in case of error
set -e

BACKEND_ENV=development \
docker-compose \
-f .deploy/shared.yml \
-f .deploy/development.yml \
config > docker-stack.yml

docker-compose -f docker-stack.yml up -d
