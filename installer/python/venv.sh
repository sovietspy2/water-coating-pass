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

# Ensure venv module is available; install it if missing
if ! python3 -c "import venv" >/dev/null 2>&1; then
    echo "[INFO] python3-venv is missing. Installing..."
    apt-get update -qq
    apt-get install -y python3-venv
fi

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