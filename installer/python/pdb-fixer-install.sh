#!/bin/bash
set -e

if ! command -v python3 >/dev/null 2>&1; then
    echo "python3 is required but not installed."
    exit 1
fi

if ! command -v pip3 >/dev/null 2>&1; then
    echo "pip3 not found. Installing python3-pip..."
    sudo apt-get update
    sudo apt-get install -y python3-pip
fi

echo "Python: $(python3 --version 2>&1)"
echo "Pip: $(pip3 --version 2>&1)"

pip3 install --upgrade pip
pip3 install pdb-fixer