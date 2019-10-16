#! /usr/bin/env bash

# Exit in case of error
set -e

DOMAIN=${DOMAIN} \
TRAEFIK_TAG=${TRAEFIK_TAG} \
STACK_NAME=${STACK_NAME} \
TAG=${TAG} \
docker-compose \
-f .deploy/shared.yml \
-f .deploy/swarm.yml \
config > docker-stack.yml

docker-auto-labels docker-stack.yml

docker stack deploy -c docker-stack.yml --with-registry-auth ${STACK_NAME}
