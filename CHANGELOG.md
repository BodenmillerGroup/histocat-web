# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.4] - 2021-06-16

- Rename `OPEN_USER_REGISTRATION` env variable according to AirLab naming scheme.
- `steinbock` default folder names updated according to v0.5.2 changes.
- Some packages updated.

## [2.1.3] - 2021-06-10

- Scatter plots optimizations. For pending issue see https://github.com/BodenmillerGroup/histocat-web/issues/17
- FIX: proper color normalization for raster images.  
- Plotly.js in now bundled instead of being served from CDN.  
- Fixed mkdocs config due to recent breaking changes in v1.2.
- Functional `steinbock` data import.

## [2.1.2] - 2021-06-07

- FIX: bug introduced in Vue 2.6.13 (https://github.com/vuejs/vue/issues/12102).
- Plotly.js updated to v2.  
- WIP: `steinbock` data import.

## [2.1.1] - 2021-06-02

- Per-channel Min-Max / z-Score normalization as a segmentation pre-processing step.

## [2.1.0] - 2021-05-27

Initial open-source release
