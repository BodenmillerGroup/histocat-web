#! /usr/bin/env bash

# Exit in case of error
set -e

source ./scripts/build.sh

docker-compose -f docker-stack.yml push
