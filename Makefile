PATH  := node_modules/.bin:$(PATH)
SHELL := /bin/bash

init:
	cd backend/app && poetry install
	cd frontend && yarn install

deploy:
	docker-compose up -d

deploy-prod:
	docker-compose -f docker/prod.yml up -d

build:
	./scripts/build.sh

build-push:
	./scripts/build-push.sh

serve:
	cd frontend && yarn serve

clean:
	sudo find . -type d -name __pycache__ -exec rm -r {} \+
	sudo find . -type d -name .pytest_cache -exec rm -r {} \+

update-frontend:
	cd frontend && ncu -u && yarn install

update-backend:
	cd backend/app && poetry update

update: update-frontend update-backend

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
