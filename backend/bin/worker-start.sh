#! /usr/bin/env bash
set -e

python3 /app/histocat/worker_pre_start.py

PROCESSES=${DRAMATIQ_PROCESSES:-8}
THREADS=${DRAMATIQ_THREADS:-8}

if [ "$BACKEND_ENV" = "development" ]
then
    echo "Starting workers in development mode"
    dramatiq --processes 1 --threads 1 --watch . histocat.worker.worker
else
    echo "Starting workers in production mode"
    dramatiq --processes $PROCESSES --threads $THREADS histocat.worker.worker
fi
