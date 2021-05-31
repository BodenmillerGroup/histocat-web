# Staging deployment

`staging` build is used to test production build of AirLab in environment similar to `production` environment, but it runs on a separate test server.
The difference is another domain name address.

Unlike production deployment, staging deployment doesn't require pre-build Docker images as it will build them during deployment process.

## Installation

Clone the histoCAT repo
```sh
git clone https://github.com/BodenmillerGroup/histocat-web.git
```

Later on you can update the repo by running:
```sh
git pull
```

Rename `.env.template` file in the root repo folder to `.env` and set all environment variables values accordingly to staging configuration.
Usually, it should be done only once as `.env` file is ignored by Git and won't be overwritten by `git pull` command.

!!! warning "Warning"
    Please check that env variable `DOMAIN` is properly set in `.env` file on staging server before deployment!

Build and deploy Docker containers on the staging server:
```sh
make deploy-staging
``` 
