#! /usr/bin/env bash

# Let the DB start
python /app/histocat/backend_pre_start.py

# Run migrations
alembic upgrade head

# Create initial data in DB
python /app/histocat/initial_data.py
