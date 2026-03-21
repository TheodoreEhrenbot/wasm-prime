#!/usr/bin/env bash
set -e
wasm-pack build --target web
elm make Main.elm --output=elm.js
