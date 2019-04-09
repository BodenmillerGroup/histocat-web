PATH  := node_modules/.bin:$(PATH)
SHELL := /bin/bash

.PHONY: clean

init:
	pip install -r backend/app/requirements/dev.txt
	cd frontend && npm install

deploy:
	docker-compose up -d

serve:
	cd frontend && npm run serve

clean:
	rm -rf build
