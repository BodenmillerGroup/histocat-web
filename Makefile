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
