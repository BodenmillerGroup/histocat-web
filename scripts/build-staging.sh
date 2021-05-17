#! /usr/bin/env bash

# Exit in case of error
set -e

if [ $(uname -s) = "Linux" ]; then
    echo "Remove __pycache__ files"
    sudo find . -type d -name __pycache__ -exec rm -r {} \+
fi

TAG=${DOCKER_TAG-latest} docker-compose -f .deploy/build-staging.yml config > docker-stack.yml

docker-compose -f docker-stack.yml build --parallel
