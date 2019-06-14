#! /usr/bin/env bash

# Exit in case of error
set -e

DOMAIN=backend \
docker-compose \
-f docker/shared.yml \
-f docker/prod.yml \
-f docker/test.yml \
config > docker-stack.yml

docker-compose -f docker-stack.yml build
docker-compose -f docker-stack.yml down -v --remove-orphans # Remove possibly previous broken stacks left hanging after an error
docker-compose -f docker-stack.yml up -d
docker-compose -f docker-stack.yml exec -T backend-tests ./tests-start.sh
docker-compose -f docker-stack.yml down -v --remove-orphans
