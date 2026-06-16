#!/usr/bin/env bash
set -Eeuo pipefail

export DEBIAN_FRONTEND=noninteractive

echo "Installing python3"
apt-get update
apt-get install -y python3