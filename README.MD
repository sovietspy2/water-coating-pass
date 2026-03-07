# PASS like coating

## Compile
- `cd src` and then `make`

## Install
- `make install`

## Pipeline scripts

The project includes two user-facing pipeline scripts in `pipeline/`:

- `pass5x1.sh`
- `pass5x5.sh`

These scripts combine:
- the `pass` executable
- the MM cleanup workflow in `mm_minim.sh`
- the water filtering step in `clean_pdb.py`

### Input for the scripts
- You need a valid `.pdb` file

### `pass5x1.sh`
Runs:
1. one PASS step with **5 water layers**
2. one MM and cleanup step

### `pass5x5.sh`
Runs:
- **5 iterations**
- each iteration does:
    - PASS with **1 water layer**
    - MM cleanup

## requirements:
- python
- linux
- C compiler

## Acknowledgments

This project reuses code from **MobyWat** (http://www.mobywat.com).]()






