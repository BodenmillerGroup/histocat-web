# Developer documentation

There are three configurations available by default: `development`, `staging` and `production`.
It is highly recommended to use `make` commands (see `Makefile` for available options) to manage all tasks.
To get a local development copy up and running please follow these steps.

## Prerequisites

Recommended and tested environment for development and deployment is Ubuntu Linux distribution. If you are going to use another OS or distribution please make changes accordingly. Make sure that the following tools are installed globally on your machine:
* [Docker](https://docs.docker.com/engine/install/ubuntu/)
* [Docker Compose](https://docs.docker.com/compose/install/)
* [Node.js LTS](https://github.com/nodesource/distributions/blob/master/README.md#debinstall)
* [Yarn Classic (v1)](https://classic.yarnpkg.com/en/docs/install#debian-stable)
* [Poetry](https://python-poetry.org/docs/#installation)

## Installation

1. Clone the repo
    ```sh
    git clone https://github.com/BodenmillerGroup/histocat-web.git
    ```
2. Install NPM packages
    ```sh
    make bootstrap
    ```
3. Deploy development Docker setup locally 
    ```sh
    make deploy-development
    ```
4. Run development version of front-end app (with hot-reloading) 
    ```sh
    make serve-frontend
    ```
5. Access locally deployed version of histoCAT-web at http://localhost:9999


## Configuration

In order to tweak the platform according to your needs, please modify environment variables configured in these files
(duplicate values across these files are possible, it should be cleaned up in the future):

1. [.env](https://github.com/BodenmillerGroup/histocat-web/blob/master/.env)
