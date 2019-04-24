PATH  := node_modules/.bin:$(PATH)
SHELL := /bin/bash

.PHONY: clean

init:
	cd backend/app && poetry install
	cd frontend && npm install

deploy:
	docker-compose up -d

docker-build:
	docker-compose build

serve:
	cd frontend && npm run serve

clean:
	rm -rf build

update-frontend:
	cd frontend && ncu -u && npm i

update-backend:
	cd backend/app && poetry update
