#!/bin/bash
set -e

# switching to src directory
cd "$(dirname "$0")/../../src"

# Compiling wdrop
echo "Compiling wdrop"
make

# Installing wdrop
echo "Installing wdrop"
make install