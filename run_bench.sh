#!/bin/sh

. $HOME/.cargo/env

echo "Running Julia code"
julia ising.jl

echo ""

echo "Running Rust code"
cargo run --release