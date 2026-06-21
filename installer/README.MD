# WDROP Installer

This installer bootstraps the local environment required to run the WDROP toolchain. It checks for required commands and Python components, runs component-specific installer scripts for anything missing, and verifies the final installation state before exiting.

## What it installs

The script verifies these components in this order:

- **Python 3** (`python3`)
- **Python 3 venv support** (usable `venv` module plus an existing local virtual environment)
- **Python module `pdbfixer`** inside the local virtual environment
- **GROMACS** (`gmx`)
- **Tinker** (`analyze`)
- **wdrop** (`wdrop`)
- **mobywat** (`mobywat`)

If any of these are missing, the installer runs the corresponding helper script:

- `gromacs/install.sh`
- `tinker/install.sh`
- `python/python3.sh`
- `python/venv.sh`
- `python/pdbfixer.sh`
- `wdrop/install.sh`
- `mobywat/install.sh`


## PATH handling

The installer ensures that the target user's `~/bin` directory exists and is available in `PATH`.

It will:

- Create `~/bin` if needed
- Create `~/.bashrc` if it does not exist
- Append this line to `~/.bashrc` if it is not already present:

```bash
export PATH="$HOME/bin:$PATH"
```

## Requirements

Before running the installer, make sure:

- You have a user who can `sudo`
- Bash is available
- The component installer scripts listed above exist in the expected directory layout
- You have permission to install system packages or software dependencies when required by the child installer scripts
- `python3-venv` support is available or installable on your system if the virtual environment needs to be created


## Example output

```text
WDROP INSTALLER

Checking required programs...

[OK] Python 3 found: /usr/bin/python3
[INSTALL] GROMACS is missing, running installer: /path/to/gromacs/install.sh

Re-checking installation state...

All required programs/modules are installed.
```

## Expected repository layout

This README assumes the installer lives in a directory structure similar to:

```text
project/
├── .venv/
├── install/
│   ├── install.sh
│   ├── gromacs/install.sh
│   ├── tinker/install.sh
│   ├── python/python3.sh
│   ├── python/venv.sh
│   ├── python/pdbfixer.sh
│   ├── wdrop/install.sh
│   └── mobywat/install.sh
```

Adjust the example paths if your top-level installer has a different filename or directory.

## Notes

- The installer is idempotent for already installed components: if a command or module is already present, it reports `[OK]` and skips the matching child installer.
- The script performs a final verification pass after all installation attempts, so partial failures are reported clearly.
- The Python module check is intentionally isolated in a subshell so the parent shell environment is not polluted by virtual environment activation.