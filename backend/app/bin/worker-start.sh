#! /usr/bin/env bash
set -e

python /app/app/worker_pre_start.py

PROCESSES=${DRAMATIQ_PROCESSES:-8}
THREADS=${DRAMATIQ_THREADS:-8}

if [ "$BACKEND_ENV" = "development" ]
then
    echo "Starting workers in development mode"
    dramatiq --processes 1 --threads 1 --watch . app.worker
else
    echo "Starting workers in production mode"
    dramatiq --processes $PROCESSES --threads $THREADS app.worker
fi
