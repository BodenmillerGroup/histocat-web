#! /usr/bin/env bash

# Exit in case of error
set -e

if [ $(uname -s) = "Linux" ]; then
    echo "Remove __pycache__ files"
    sudo find . -type d -name __pycache__ -exec rm -r {} \+
fi

BACKEND_ENV=production \
docker-compose -f .deploy/build.yml config > docker-stack.yml

docker-compose -f docker-stack.yml build --parallel
