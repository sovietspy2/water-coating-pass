#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="${SCRIPT_DIR}/../../.venv"
TEST_VENV="/tmp/test-venv-check.$$"

echo "[INFO] script dir: ${SCRIPT_DIR}"
echo "[INFO] venv dir:   ${VENV_DIR}"

if ! command -v python3 >/dev/null 2>&1; then
    echo "[ERROR] python3 is required but not installed."
    exit 1
fi

PYVER="$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"
echo "[INFO] Detected Python version: ${PYVER}"

rm -rf "$TEST_VENV"
if ! python3 -m venv "$TEST_VENV" >/dev/null 2>&1; then
    echo "[INFO] venv support missing or broken, installing matching package..."
    sudo apt update -qq
    sudo apt install -y "python${PYVER}-venv" python3-pip
    rm -rf "$TEST_VENV"
    python3 -m venv "$TEST_VENV" >/dev/null 2>&1
fi
rm -rf "$TEST_VENV"

if [ -f "$VENV_DIR/bin/activate" ]; then
    echo "[OK] Virtual environment already exists at $VENV_DIR"
else
    echo "[INFO] Creating virtual environment at $VENV_DIR..."
    rm -rf "$VENV_DIR"
    python3 -m venv "$VENV_DIR"
    echo "[OK] Virtual environment created."
fi

(
    source "$VENV_DIR/bin/activate"
    echo "[INFO] Python: $(python --version)"
    echo "[INFO] Pip:    $(python -m pip --version)"
    python -m pip install --upgrade pip --quiet
)