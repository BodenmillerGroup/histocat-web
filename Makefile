PATH  := node_modules/.bin:$(PATH)
SHELL := /bin/bash

bootstrap:
	cd backend && JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64 poetry install --extras "backend worker"
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

build-push-worker-base-image:
	./scripts/build-worker-base-image.sh

clean:
	sudo find . -type d -name __pycache__ -exec rm -r {} \+
	find . -type d -name .pytest_cache -exec rm -r {} \+
	find . -type d -name .mypy_cache -exec rm -r {} \+
	rm docker-stack.yml

serve-frontend:
	cd frontend && yarn run serve

update-frontend:
	cd frontend && ncu -u && yarn install

update-backend:
	cd backend && poetry update

mypy:
	cd backend && mypy histocat

pyright:
	cd backend && pyright

isort:
	cd backend && isort --multi-line=3 --trailing-comma --force-grid-wrap=0 --combine-as --line-width 88 histocat

vulture:
	cd backend && vulture histocat --min-confidence 70

black:
	cd backend && black histocat

autoflake:
	cd backend && autoflake --remove-all-unused-imports --recursive --remove-unused-variables --in-place app --exclude=__init__.py
