#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="${SCRIPT_DIR}/../../.venv"

echo "[INFO] script dir: ${SCRIPT_DIR}"
echo "[INFO] venv dir:   ${VENV_DIR}"

if ! command -v python3 >/dev/null 2>&1; then
    echo "[ERROR] python3 is required but not installed."
    exit 1
fi

PYVER="$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"

if ! python3 -m venv --without-pip /tmp/test-venv-check.$$ >/dev/null 2>&1; then
    echo "[INFO] venv support looks broken. Installing matching package..."
    sudo apt update -qq
    sudo apt install -y "python${PYVER}-venv"
fi

rm -rf /tmp/test-venv-check.$$

if [ -d "$VENV_DIR" ] && [ -f "$VENV_DIR/bin/activate" ]; then
    echo "[OK] Virtual environment already exists at $VENV_DIR"
else
    echo "[INFO] Creating virtual environment at $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
    echo "[OK] Virtual environment created."
fi

# Upgrade pip inside the venv using a subshell so activation doesn't leak
(
    # shellcheck source=/dev/null
    source "$VENV_DIR/bin/activate"
    echo "[INFO] Python: $(python3 --version)"
    echo "[INFO] Pip:    $(python3 -m pip --version)"
    python3 -m pip install --upgrade pip --quiet
)