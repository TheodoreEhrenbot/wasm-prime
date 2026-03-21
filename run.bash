#!/usr/bin/env bash
set -e
./build.bash
python3 -m http.server
