#!/usr/bin/env bash
set -Eeuo pipefail

export DEBIAN_FRONTEND=noninteractive

PYTHON_CMD=""

check_python_is_v3() {
    "$1" - <<'PY' >/dev/null 2>&1
import sys
raise SystemExit(0 if sys.version_info[0] == 3 else 1)
PY
}

detect_python() {
    if command -v python3 >/dev/null 2>&1; then
        PYTHON_CMD="python3"
        return 0
    fi

    if command -v python >/dev/null 2>&1 && check_python_is_v3 python; then
        PYTHON_CMD="python"
        return 0
    fi

    return 1
}

install_python_stack() {
    echo "Python 3 not found. Installing python3, pip, and venv..."
    apt-get update
    apt-get install -y python3 python3-pip python3-venv
}

ensure_pip() {
    if "$PYTHON_CMD" -m pip --version >/dev/null 2>&1; then
        return 0
    fi

    echo "pip for Python 3 not found. Trying ensurepip..."
    if "$PYTHON_CMD" -m ensurepip --upgrade >/dev/null 2>&1; then
        if "$PYTHON_CMD" -m pip --version >/dev/null 2>&1; then
            return 0
        fi
    fi

    echo "ensurepip did not provide pip. Installing python3-pip via apt..."
    apt-get update
    apt-get install -y python3-pip python3-venv

    "$PYTHON_CMD" -m pip --version >/dev/null 2>&1
}

pip_upgrade_with_fallback() {
    local package="$1"

    if ! "$PYTHON_CMD" -m pip install --upgrade "$package"; then
        echo "Retrying ${package} install with PIP_BREAK_SYSTEM_PACKAGES=1..."
        PIP_BREAK_SYSTEM_PACKAGES=1 "$PYTHON_CMD" -m pip install --upgrade "$package"
    fi
}

if ! detect_python; then
    install_python_stack
    detect_python
fi

ensure_pip

echo "Using interpreter: $("$PYTHON_CMD" -c 'import sys; print(sys.executable)')"
echo "Detected version: $("$PYTHON_CMD" --version 2>&1)"
echo "Detected pip: $("$PYTHON_CMD" -m pip --version 2>&1)"

pip_upgrade_with_fallback pip
pip_upgrade_with_fallback pdbfixer

echo "python installation completed successfully."