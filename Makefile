PATH  := node_modules/.bin:$(PATH)
SHELL := /bin/bash

.PHONY: clean

init:
	pip install -r backend/app/requirements/dev.txt
	cd frontend && npm install

deploy:
	docker-compose up -d

docker-build:
	docker-compose build

serve:
	cd frontend && npm run serve

clean:
	rm -rf build

upgrade-frontend:
	cd frontend && ncu -u && npm i

upgrade-backend-dev:
	pip-compile --upgrade backend/app/requirements/dev.in

upgrade-backend-prod:
	pip-compile --upgrade backend/app/requirements/prod.in

sync:
	pip-sync backend/app/requirements/dev.txt
