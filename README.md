# WebGestalt Rust

[![Rust](https://github.com/bzhanglab/webgestalt_rust/actions/workflows/rust.yml/badge.svg?branch=master)](https://github.com/bzhanglab/webgestalt_rust/actions/workflows/rust.yml)

Rust implementation of [WebGestaltR](https://github.com/bzhanglab/webgestaltr).

## Notes

This CLI is focused purely on computation. **It does not provide GMT files or HTML reports**. The output of this tool is JSON files containing the results. For a more feature-complete tool, see the original [WebGestaltR](https://bzhanglab.github.io/WebGestaltR/) tool.

## Install

Requires the [rust toolchain](https://rustup.rs/). Then run the following command in your terminal:

```shell
cargo install webgestalt
```

## CLI

For help with CLI, run

```shell
webgestalt --help
```

Example of running over-representation analysis using `kegg.gmt`, with an interesting list at `int.txt` and a reference of `ref.txt`. Outputs JSON file at `output.json`

```shell
webgestalt ora -g kegg.gmt -i int.txt -r ref.txt -o output.json
```