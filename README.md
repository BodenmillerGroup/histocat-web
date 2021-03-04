<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->


<!-- ABOUT THE PROJECT -->
# histoCAT-web

[![histocat][product-screenshot]](https://github.com/BodenmillerGroup/histocat-web)

Histology Topography Cytometry Analysis Toolbox (histoCAT) is a web platform to visualize and analyse multiplexed image cytometry data interactively.

<!-- GETTING STARTED -->
## Getting Started

There are three configurations available by default: `development`, `staging` and `production`.
It is highly recommended to use `make` commands (see `Makefile` for available options) to manage all tasks.
To get a local development copy up and running please follow these steps.

### Prerequisites

Recommended and tested environment for development and deployment is Ubuntu Linux distribution. If you are going to use another OS or distribution please make changes accordingly. Make sure that the following tools are installed globally on your machine:
* [Docker](https://docs.docker.com/engine/install/ubuntu/)
* [Docker Compose](https://docs.docker.com/compose/install/)
* [Node.js LTS](https://github.com/nodesource/distributions/blob/master/README.md#debinstall)
* [Yarn Classic (v1)](https://classic.yarnpkg.com/en/docs/install#debian-stable)
* [Poetry](https://python-poetry.org/docs/#installation)

### Installation

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


### Configuration

In order to tweak the platform according to your needs, please modify environment variables configured in these files
(duplicate values across these files are possible, it should be cleaned up in the future):
1. [.env](.env)
2. [packages/frontend/.env](packages/frontend/.env)


<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/BodenmillerGroup/histocat-web/issues) for a list of proposed features (and known issues).



<!-- CONTRIBUTING -->
## Contributing

Any contributions you make are greatly appreciated.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request


<!-- CONTACT -->
## Contact

Anton Rau - anton.rau@uzh.ch

Project Link: [https://github.com/BodenmillerGroup/histocat-web](https://github.com/BodenmillerGroup/histocat-web)



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[product-screenshot]: frontend/public/img/histoCAT.png
