#!/bin/bash
set -e
echo "WDROP INSTALLER"

# functions
run_installer() {
    local name="$1"
    local script_path="$2"

    echo "[INSTALL] $name is missing, running installer: $script_path"

    if [ ! -x "$script_path" ]; then
        echo "[ERROR] Installer script is not executable or does not exist: $script_path"
        exit 1
    fi

    "$script_path"
    installed_any=1
}

check_cmd() {
    local label="$1"
    local cmd="$2"

    if command -v "$cmd" >/dev/null 2>&1; then
        echo "[OK] $label found: $(command -v "$cmd")"
    else
        echo "[MISSING] $label ($cmd) not found"
        missing=1
    fi
}

check_python_module() {
    local module="$1"

    if command -v python3 >/dev/null 2>&1; then
        if python3 -c "import ${module}" >/dev/null 2>&1; then
            echo "[OK] Python module found: $module"
        else
            echo "[MISSING] Python module not found: $module"
            missing=1
        fi
    else
        echo "[MISSING] python3 not found, so module check skipped: $module"
        missing=1
    fi
}

# Ading ~/bin to path
mkdir -p "$HOME/bin"
grep -qxF 'export PATH="$HOME/bin:$PATH"' "$HOME/.bashrc" || echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
export PATH="$HOME/bin:$PATH"


# Install script locations
INSTALL_GROMACS_SCRIPT="./gromacs/install.sh"
INSTALL_TINKER_SCRIPT="./tinker/install.sh"
INSTALL_PYTHON3_SCRIPT="/python/install.sh"
INSTALL_PDBFIXER_SCRIPT="./python/pdb-fixer-install.sh"
INSTALL_WDROP_SCRIPT="./wdrop/install.sh"
INSTALL_MOBYWAT_SCRIPT="./mobywat/install.sh"

missing_any=0
installed_any=0

echo "Checking required programs..."
echo

if ! command -v gmx >/dev/null 2>&1; then
    missing_any=1
    run_installer "GROMACS" "$INSTALL_GROMACS_SCRIPT"
else
    echo "[OK] GROMACS found"
fi

if ! command -v analyze >/dev/null 2>&1; then
    missing_any=1
    run_installer "Tinker" "$INSTALL_TINKER_SCRIPT"
else
    echo "[OK] Tinker found"
fi

if ! command -v python3 >/dev/null 2>&1; then
    missing_any=1
    run_installer "Python 3" "$INSTALL_PYTHON3_SCRIPT"
else
    echo "[OK] Python 3 found"
fi

if command -v python3 >/dev/null 2>&1; then
    if ! python3 -c "import pdbfixer" >/dev/null 2>&1; then
        missing_any=1
        run_installer "Python module pdbfixer" "$INSTALL_PDBFIXER_SCRIPT"
    else
        echo "[OK] Python module pdbfixer found"
    fi
else
    echo "[WARN] Skipping pdbfixer check because python3 is not available yet"
fi

if ! command -v wdrop >/dev/null 2>&1; then
    missing_any=1
    run_installer "wdrop" "$INSTALL_WDROP_SCRIPT"
else
    echo "[OK] wdrop found"
fi

if ! command -v mobywat >/dev/null 2>&1; then
    missing_any=1
    run_installer "mobywat" "$INSTALL_MOBYWAT_SCRIPT"
else
    echo "[OK] mobywat found"
fi

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

if ! command -v gmx >/dev/null 2>&1; then
    missing_items+=("GROMACS (gmx)")
fi

if ! command -v analyze >/dev/null 2>&1; then
    missing_items+=("Tinker")
fi

if ! command -v python3 >/dev/null 2>&1; then
    missing_items+=("Python 3")
fi

if command -v python3 >/dev/null 2>&1; then
    if ! python3 -c "import pdbfixer" >/dev/null 2>&1; then
        missing_items+=("Python module pdbfixer")
    fi
else
    missing_items+=("Python module pdbfixer (python3 unavailable)")
fi

if ! command -v wdrop >/dev/null 2>&1; then
    missing_items+=("wdrop")
fi

if ! command -v mobywat >/dev/null 2>&1; then
    missing_items+=("wdrop")
fi

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