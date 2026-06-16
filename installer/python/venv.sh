#!/bin/bash
set -e

VENV_DIR="../.venv" # project root

if ! command -v python3 >/dev/null 2>&1; then
    echo "python3 is required but not installed."
    exit 1
fi

if ! python3 -m venv --help >/dev/null 2>&1; then
    echo "python3 venv support is missing. Installing python3-venv..."
    sudo apt-get update
    sudo apt-get install -y python3-venv
fi

if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment in $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
fi

# Activate venv
source "$VENV_DIR/bin/activate"

echo "Virtual environment activated."
echo "Python: $(python --version 2>&1)"
echo "Pip: $(python -m pip --version 2>&1)"