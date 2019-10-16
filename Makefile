PATH  := node_modules/.bin:$(PATH)
SHELL := /bin/bash

init:
	cd backend/app && poetry install --extras backend
	cd frontend && yarn install

deploy-development:
	./scripts/deploy-development.sh

deploy-production:
	./scripts/deploy-production.sh

deploy-staging:
	./scripts/deploy-staging.sh

build:
	./scripts/build.sh

build-push:
	./scripts/build-push.sh

serve:
	cd frontend && yarn serve

clean:
	sudo find . -type d -name __pycache__ -exec rm -r {} \+
	find . -type d -name .pytest_cache -exec rm -r {} \+
	find . -type d -name .mypy_cache -exec rm -r {} \+
	rm docker-stack.yml

update-frontend:
	cd frontend && ncu -u && yarn install

update-backend:
	cd backend/app && poetry update

mypy:
	cd backend/app && mypy app

pyright:
	cd backend/app && pyright

isort:
	cd backend/app && isort --multi-line=3 --trailing-comma --force-grid-wrap=0 --combine-as --line-width 88 --recursive --apply app

vulture:
	cd backend/app && vulture app --min-confidence 70

black:
	cd backend/app && black app

autoflake:
	cd backend/app && autoflake --remove-all-unused-imports --recursive --remove-unused-variables --in-place app --exclude=__init__.py
