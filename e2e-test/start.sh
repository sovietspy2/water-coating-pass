#!/usr/bin/env bash
set -Eeuo pipefail

IMAGE_NAME="wdrop"
HOST_OUTPUT_DIR="${HOST_OUTPUT_DIR:-$PWD/test_output}"
CONTAINER_OUTPUT_DIR="/data/test_runs"

mkdir -p "$HOST_OUTPUT_DIR"

docker build -t "$IMAGE_NAME" .

docker run --rm -it \
-v "${HOST_OUTPUT_DIR}:${CONTAINER_OUTPUT_DIR}" \
"$IMAGE_NAME" \
/bin/bash -lc "cd /opt/water-coating-pass/testing && ./test.sh -p ${CONTAINER_OUTPUT_DIR} $*"