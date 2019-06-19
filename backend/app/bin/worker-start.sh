#! /usr/bin/env bash
set -e

python /app/app/worker_pre_start.py

if [ "$BACKEND_ENV" = "development" ]
then
    echo "Starting workers in development mode"
    dramatiq --processes 1 --threads 1 --watch . app.worker
else
    echo "Starting workers in production mode"
    dramatiq --processes 1 --threads 8 app.worker
fi
