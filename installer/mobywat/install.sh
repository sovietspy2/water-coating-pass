#!/bin/bash
set -e

URL="http://mobywat.com/files/download/ver_11/mobywat.tgz"
TMPDIR="$(mktemp -d /tmp/mobywat.XXXXXX)"
ARCHIVE="$TMPDIR/mobywat.tgz"

cleanup() {
    rm -rf "$TMPDIR"
}

trap cleanup EXIT INT TERM

echo "Downloading mobywat source..."
if command -v curl >/dev/null 2>&1; then
    curl -L "$URL" -o "$ARCHIVE"
elif command -v wget >/dev/null 2>&1; then
    wget -O "$ARCHIVE" "$URL"
else
    echo "Error: neither curl nor wget is installed."
    exit 1
fi

echo "Extracting archive to $TMPDIR ..."
tar -xzf "$ARCHIVE" -C "$TMPDIR"

echo "Locating src directory..."
SRC_DIR="$(find "$TMPDIR" -type d -name src | head -n 1)"

if [ -z "$SRC_DIR" ]; then
    echo "Error: could not find src directory after extraction."
    exit 1
fi

echo "Entering $SRC_DIR"
cd "$SRC_DIR"

echo "Running make..."
make

echo "Running make install..."
make install

echo "Installation completed."