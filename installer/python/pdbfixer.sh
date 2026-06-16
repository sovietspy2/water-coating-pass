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

if ! python3 -c "import venv" >/dev/null 2>&1; then
    echo "[ERROR] python3-venv is required but not available."
    echo "[ERROR] Run the venv installer first: installer/python/venv.sh"
    exit 1
fi

# Create the venv if it doesn't exist yet (defensive: venv.sh should have run first)
if [ ! -d "$VENV_DIR" ] || [ ! -f "$VENV_DIR/bin/activate" ]; then
    echo "[INFO] Venv not found, creating it at $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
fi

# Install pdbfixer inside the venv using a subshell
(
    # shellcheck source=/dev/null
    source "$VENV_DIR/bin/activate"
    echo "[INFO] Python: $(python3 --version)"
    echo "[INFO] Pip:    $(python3 -m pip --version)"
    python3 -m pip install --upgrade pip --quiet
    python3 -m pip install pdbfixer
    echo "[OK] pdbfixer installed: $(python3 -c 'import pdbfixer; print(pdbfixer.__file__)')"
)