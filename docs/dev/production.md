# Production deployment

## Preparing production Docker images

Production deployment requires pre-build production Docker images to be available.
Images names are defined by the following env variables:

- DOCKER_IMAGE_BACKEND
- DOCKER_IMAGE_WORKER
- DOCKER_IMAGE_FRONTEND

!!! important "Important"
    It is strongly recommended to use tagged Docker images in production deployment scenario to enforce stability! `DOCKER_TAG` env variable in `.env` file defines Docker image tag used during build and deployment. 

In order to build these images and then push them to Docker image registry (e.g. Docker Hub), run the following command in your development environment:
```sh
make build-push
```

This operation should be done each time you are planning to release new production version.

## Installation

As soon as Docker images are build and published to a registry, one can deploy production version of histoCAT by running the following commands on a production server:

Clone the histoCAT repo
```sh
git clone https://github.com/BodenmillerGroup/histocat-web.git
```

Later on you can update the repo by running:
```sh
git pull
```

Rename `.env.template` file in the root repo folder to `.env` and set all environment variables values accordingly to production configuration.
Usually it should be done only once as `.env` file is ignored by Git and won't be overwritten by `git pull` command.

!!! warning "Warning"
    Please check that env variable `DOMAIN` is properly set in `.env` file on production server before deployment!

Deploy production Docker containers 
```sh
make deploy-production
```

!!! info "Info"
    Production deployment doesn't support debugging and hot-reloading because it is compiled to minimize bundle size and optimize performance. 
