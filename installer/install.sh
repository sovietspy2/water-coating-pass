#!/usr/bin/env bash
set -Eeuo pipefail

echo "WDROP INSTALLER"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

TARGET_USER="${SUDO_USER:-${USER:-root}}"
TARGET_HOME="${HOME}"

if [ -n "${SUDO_USER:-}" ] && command -v getent >/dev/null 2>&1; then
    TARGET_HOME="$(getent passwd "$SUDO_USER" | cut -d: -f6)"
fi

HOME_BIN="${TARGET_HOME}/bin"
BASHRC_FILE="${TARGET_HOME}/.bashrc"
PATH_LINE='export PATH="$HOME/bin:$PATH"'

missing_any=0
installed_any=0

ensure_root() {
    if [ "${EUID:-$(id -u)}" -eq 0 ]; then
        return 0
    fi

    if command -v sudo >/dev/null 2>&1; then
        echo "[INFO] Re-running with sudo..."
        exec sudo -E bash "$0" "$@"
    fi

    echo "[ERROR] Root privileges are required, but sudo is not available."
    echo "[ERROR] Run this script as root."
    echo "[ERROR] In Docker, enter the container as root instead of relying on sudo."
    exit 1
}

add_home_bin_to_path() {
    mkdir -p "$HOME_BIN"

    if [ ! -f "$BASHRC_FILE" ]; then
        touch "$BASHRC_FILE"
    fi

    grep -qxF "$PATH_LINE" "$BASHRC_FILE" || printf '%s\n' "$PATH_LINE" >> "$BASHRC_FILE"

    case ":$PATH:" in
        *":$HOME_BIN:"*) ;;
        *) export PATH="$HOME_BIN:$PATH" ;;
    esac
}

run_installer() {
    local name="$1"
    local script_path="$2"

    echo "[INSTALL] $name is missing, running installer: $script_path"

    if [ ! -f "$script_path" ]; then
        echo "[ERROR] Installer script does not exist: $script_path"
        exit 1
    fi

    bash "$script_path"
    installed_any=1
}

has_cmd() {
    command -v "$1" >/dev/null 2>&1
}

has_python3() {
    has_cmd python3
}

has_python_venv() {
    has_python3 || return 1
    python3 -m venv --help >/dev/null 2>&1
}

has_python_module() {
    local module="$1"
    has_python3 || return 1
    python3 -c "import ${module}" >/dev/null 2>&1
}

ensure_component_cmd() {
    local label="$1"
    local cmd="$2"
    local installer="$3"

    if has_cmd "$cmd"; then
        echo "[OK] $label found: $(command -v "$cmd")"
    else
        missing_any=1
        run_installer "$label" "$installer"
    fi
}

ensure_component_venv() {
    local label="$1"
    local installer="$2"

    if has_python_venv; then
        echo "[OK] $label found"
    else
        missing_any=1
        run_installer "$label" "$installer"
    fi
}

ensure_component_module() {
    local label="$1"
    local module="$2"
    local installer="$3"

    if has_python_module "$module"; then
        echo "[OK] $label found"
    else
        missing_any=1
        run_installer "$label" "$installer"
    fi
}

collect_missing_cmd() {
    local label="$1"
    local cmd="$2"

    if ! has_cmd "$cmd"; then
        missing_items+=("$label")
    fi
}

collect_missing_venv() {
    local label="$1"

    if ! has_python_venv; then
        missing_items+=("$label")
    fi
}

collect_missing_module() {
    local label="$1"
    local module="$2"

    if ! has_python3; then
        missing_items+=("${label} (python3 unavailable)")
    elif ! has_python_module "$module"; then
        missing_items+=("$label")
    fi
}

INSTALL_GROMACS_SCRIPT="${SCRIPT_DIR}/gromacs/install.sh"
INSTALL_TINKER_SCRIPT="${SCRIPT_DIR}/tinker/install.sh"
INSTALL_PYTHON3_SCRIPT="${SCRIPT_DIR}/python/python3.sh"
INSTALL_PYTHON3_VENV="${SCRIPT_DIR}/python/venv.sh"
INSTALL_PDBFIXER_SCRIPT="${SCRIPT_DIR}/python/pdb-fixer.sh"
INSTALL_WDROP_SCRIPT="${SCRIPT_DIR}/wdrop/install.sh"
INSTALL_MOBYWAT_SCRIPT="${SCRIPT_DIR}/mobywat/install.sh"

ensure_root "$@"
add_home_bin_to_path

echo
echo "Checking required programs..."
echo

ensure_component_cmd "GROMACS" "gmx" "$INSTALL_GROMACS_SCRIPT"
ensure_component_cmd "Tinker" "analyze" "$INSTALL_TINKER_SCRIPT"
ensure_component_cmd "Python 3" "python3" "$INSTALL_PYTHON3_SCRIPT"
ensure_component_venv "Python 3 venv" "$INSTALL_PYTHON3_VENV"
ensure_component_module "Python module pdbfixer" "pdbfixer" "$INSTALL_PDBFIXER_SCRIPT"
ensure_component_cmd "wdrop" "wdrop" "$INSTALL_WDROP_SCRIPT"
ensure_component_cmd "mobywat" "mobywat" "$INSTALL_MOBYWAT_SCRIPT"

echo
if [ "$missing_any" -eq 0 ]; then
    echo "All required programs/modules are already installed."
else
    echo "One or more missing components were detected and installer scripts were executed."
fi

echo
echo "Re-checking installation state..."
echo

missing_items=()

collect_missing_cmd "GROMACS (gmx)" "gmx"
collect_missing_cmd "Tinker (analyze)" "analyze"
collect_missing_cmd "Python 3" "python3"
collect_missing_venv "Python 3 venv"
collect_missing_module "Python module pdbfixer" "pdbfixer"
collect_missing_cmd "wdrop" "wdrop"
collect_missing_cmd "mobywat" "mobywat"

if [ "${#missing_items[@]}" -eq 0 ]; then
    echo "All required programs/modules are installed."
    exit 0
else
    echo "Some required programs/modules are still missing after installation attempts:"
    for item in "${missing_items[@]}"; do
        echo "  - $item"
    done
    exit 1
fi