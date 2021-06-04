# Development deployment

Local development deployment supports debugging and hot-reloading both on back-end and front-end.
To get a local development instance up and running please follow the following steps.

## Prerequisite

Recommended and tested environment for development and deployment is Ubuntu Linux distribution.
If you are going to use another OS or distribution please make changes accordingly.
Make sure that the following tools installed globally on your machine:

* [Git](https://git-scm.com/)
* [GNU Make](https://www.gnu.org/software/make/)
* [Docker](https://docs.docker.com/engine/install/ubuntu/)
* [Docker Compose](https://docs.docker.com/compose/install/)
* [Node.js LTS](https://github.com/nodesource/distributions/blob/master/README.md#debinstall)
* [Yarn Classic (v1)](https://classic.yarnpkg.com/en/docs/install#debian-stable)
* [Poetry](https://python-poetry.org/docs/#installation)

!!! info "Info"
    For local development it is strongly advised to use [pyenv](https://github.com/pyenv/pyenv) tool to manage an exact Python version that is used by back-end Docker containers.
    As time of the writing, it is Python 3.8.10. It is also suggested to use Python virtual environments (see [Creation of virtual environments](https://docs.python.org/3.8/library/venv.html)).

## Installation

Clone the histoCAT repo
```sh
git clone https://github.com/BodenmillerGroup/histocat-web.git
```

Install all necessary PyPI and NPM packages on your local machine
```sh
make bootstrap
```

Deploy development Docker containers locally 
```sh
make deploy-development
```

One can access locally deployed version of AirLab at [http://localhost](http://localhost).

Run development version of front-end app (with automatic hot-reloading).
```sh
make serve-frontend
```

One can access locally deployed version of AirLab with hot-reloading and debugging capabilities at [http://localhost:9999](http://localhost:9999).
Please keep in mind that this functionality is only supported on local machines with `DOMAIN` env variable set to `localhost`!
