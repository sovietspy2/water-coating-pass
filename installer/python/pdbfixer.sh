#!/bin/bash
set -e

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="${SCRIPT_DIR}/../../.venv"

echo "script dir: ${SCRIPT_DIR}"
echo "script dir: ${VENV_DIR}"

if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment in $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
fi

source "$VENV_DIR/bin/activate"

if ! command -v python3 >/dev/null 2>&1; then
    echo "python3 is required but not installed."
    exit 1
fi

if ! command -v pip3 >/dev/null 2>&1; then
    echo "pip3 not found. Installing python3-pip..."
    apt-get update
    apt-get install -y python3-pip
fi

echo "Python: $(python3 --version 2>&1)"
echo "Pip: $(pip3 --version 2>&1)"

pip3 install --upgrade pip
pip3 install pdbfixer