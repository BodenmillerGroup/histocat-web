#! /usr/bin/env bash
set -e

python /app/histocat/tests_pre_start.py

pytest $* /app/histocat/tests/
