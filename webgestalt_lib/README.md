# webgestalt_lib

Library for performing different types of enrichment analyses. Serves as the underlying core for [WebGestaltPy](https://github.com/bzhanglab/webgestaltpy), [WebGestaltR](https://github.com/bzhanglab/webgestaltr), and the pure-Rust CLI in this repository.

## Methods

Supported methods include:

- Over-representation analysis (ORA)
- Gene Set Enrichment Analysis (GSEA)

## Installation

To use webgestalt_lib in your Rust project, add the following line to your `Cargo.toml`.

```toml
webgestalt_lib = {git = "https://github.com/bzhanglab/webgestalt_rust.git"}
```

If you are just interested in running an analysis, rather than develop new tools, please use on of the packages mentioned at the beginning of the README.

## Development Priorities

1. Fast and correct implementations of enrichment methods
2. Full compatibility with the WebGestaltR package
   - The R package provides the most reporting functionality, and project was initially created to only assist the R package with the computation aspects
3. Fast compilation times
   - Every package install has to build the library from scratch, so the lower number of dependencies, the better

This crate does not provide any data formatting, or charts to display the results of the analysis. This work has already been done by the [R package](https://github.com/bzhanglab/webgestaltr), and a limited implementation is provided by the Rust CLI. The focus for this library is purely computational.
