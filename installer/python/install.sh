#!/bin/bash
set -e

export DEBIAN_FRONTEND=noninteractive

if [ "${EUID}" -eq 0 ]; then
    SUDO=""
else
    SUDO="sudo"
fi

PYTHON_CMD=""

check_python_is_v3() {
    "$1" - <<'PY' >/dev/null 2>&1
import sys
raise SystemExit(0 if sys.version_info[0] == 3 else 1)
PY
}

if command -v python3 >/dev/null 2>&1; then
    PYTHON_CMD="python3"
elif command -v python >/dev/null 2>&1 && check_python_is_v3 python; then
    PYTHON_CMD="python"
fi

if [ -z "$PYTHON_CMD" ]; then
    echo "Python 3 not found. Installing python3, pip, and venv..."
    $SUDO apt-get update
    $SUDO apt-get install -y python3 python3-pip python3-venv
    PYTHON_CMD="python3"
fi

ensure_pip() {
    if "$PYTHON_CMD" -m pip --version >/dev/null 2>&1; then
        return 0
    fi

    echo "pip for Python 3 not found. Trying ensurepip..."
    if "$PYTHON_CMD" -m ensurepip --upgrade >/dev/null 2>&1; then
        "$PYTHON_CMD" -m pip --version >/dev/null 2>&1 && return 0
    fi

    echo "ensurepip did not provide pip. Installing python3-pip via apt..."
    $SUDO apt-get update
    $SUDO apt-get install -y python3-pip python3-venv

    "$PYTHON_CMD" -m pip --version >/dev/null 2>&1
}

ensure_pip

echo "Using interpreter: $("$PYTHON_CMD" -c 'import sys; print(sys.executable)')"
echo "Detected version: $("$PYTHON_CMD" --version 2>&1)"
echo "Detected pip: $("$PYTHON_CMD" -m pip --version 2>&1)"

if ! "$PYTHON_CMD" -m pip install --upgrade pip; then
    echo "Retrying pip upgrade with --break-system-packages..."
    "$PYTHON_CMD" -m pip install --break-system-packages --upgrade pip
fi

if ! "$PYTHON_CMD" -m pip install --upgrade pdbfixer; then
    echo "Retrying pdbfixer install with --break-system-packages..."
    "$PYTHON_CMD" -m pip install --break-system-packages --upgrade pdbfixer
fi

echo "python & pdb-fixer installation completed successfully."