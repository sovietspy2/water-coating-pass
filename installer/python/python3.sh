#!/usr/bin/env bash
set -Eeuo pipefail

export DEBIAN_FRONTEND=noninteractive

echo "[INSTALL] Installing python3..."
apt-get update -qq
apt-get install -y python3 python3-venv python3-pip

echo "[OK] Python 3: $(python3 --version)"