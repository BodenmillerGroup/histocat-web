#! /usr/bin/env bash

# Exit in case of error
set -e

DOMAIN=${DOMAIN} \
TRAEFIK_TAG=${TRAEFIK_TAG} \
STACK_NAME=${STACK_NAME} \
TAG=${TAG} \
docker-compose \
-f docker/shared.admin.yml \
-f docker/shared.base-images.yml \
-f docker/shared.depends.yml \
-f docker/shared.env.yml \
-f docker/deploy.command.yml \
-f docker/deploy.images.yml \
-f docker/deploy.labels.yml \
-f docker/deploy.networks.yml \
-f docker/deploy.volumes-placement.yml \
config > docker-stack.yml

docker-auto-labels docker-stack.yml

docker stack deploy -c docker-stack.yml --with-registry-auth ${STACK_NAME}
