#! /usr/bin/env bash

# Exit in case of error
set -e

TAG=${DOCKER_TAG-latest} source ./scripts/build.sh

docker-compose -f docker-stack.yml push
