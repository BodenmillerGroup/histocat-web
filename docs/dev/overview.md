# Getting started

histoCAT application is split into three parts:

- **backend**: main API gateway. Deployed in a separate container.
- **worker**: background worker that handles long-running tasks, like data processing, image segmentation, sending emails, etc. Deployed in a separate container.
- **frontend**: front-end Vue application, served as a static content by nginx server.

There are other services deployed together with histoCAT application:

- Postgres database
- Redis (cache/pub-sub)
- RabbitMQ (job queue)
- Traefik (load-balancer)
- pgAdmin (Postgres admin app)
- RedisInsight (Redis admin app)
- Portainer (Docker admin app)

!!! info "Info"
    Services configuration is defined in `.deploy/shared.yml` file. Please also check `.deploy/development.yml`, `.deploy/production.yml` or `.deploy/staging.yml` depending on your deployment scenario.

Three deployment configurations are available by default: `development`, `staging` and `production`.
It is highly recommended to use `make` commands (see `Makefile` for available options) to manage all tasks.

By default, there is a [.env.template](https://github.com/BodenmillerGroup/histocat-web/blob/master/.env.template) file in the repo with default values.
Please rename it to `.env` and set all environment variables values according to your needs and server configuration.

!!! warning "Warning"
    `.env` file is ignored by Git in order to avoid putting sensitive information to Git repository!

!!! important "Important"
    Main setting you should care about is `DOMAIN` environment variable. In case of local development (with hot-reloading dev proxy server) it should be set to `localhost`.
