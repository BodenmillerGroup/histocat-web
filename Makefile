PATH  := node_modules/.bin:$(PATH)
SHELL := /bin/bash

.PHONY: clean

init:
	cd backend/app && poetry install
	cd frontend && yarn install

deploy:
	docker-compose up -d

docker-build:
	docker-compose build

serve:
	cd frontend && yarn serve

clean:
	rm -rf build

update-frontend:
	cd frontend && ncu -u && yarn install

update-backend:
	cd backend/app && poetry update

mypy:
	cd backend/app && mypy app

isort:
	cd backend/app && isort --multi-line=3 --trailing-comma --force-grid-wrap=0 --combine-as --line-width 88 --recursive --apply app

vulture:
	cd backend/app && vulture app --min-confidence 70

black:
	cd backend/app && black app

autoflake:
	cd backend/app && autoflake --remove-all-unused-imports --recursive --remove-unused-variables --in-place app --exclude=__init__.py
