#!/usr/bin/env bash
set -euo pipefail

GROMACS_VERSION="2026.1"
PREFIX="$HOME/opt/gromacs"
SRC_DIR="$HOME/src"
BUILD_ROOT="$SRC_DIR/gromacs-build"
TARBALL="gromacs-${GROMACS_VERSION}.tar.gz"
URL="https://ftp.gromacs.org/gromacs/${TARBALL}"
BIN_DIR="$HOME/bin"

echo "==> Installing build dependencies"
sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    git \
    wget \
    ca-certificates \
    pkg-config

echo "==> Creating directories"
mkdir -p "$SRC_DIR" "$BUILD_ROOT" "$PREFIX" "$BIN_DIR"

cd "$SRC_DIR"

echo "==> Downloading GROMACS ${GROMACS_VERSION}"
wget -O "$TARBALL" "$URL"

echo "==> Extracting source"
tar xfz "$TARBALL"

cd "gromacs-${GROMACS_VERSION}"
mkdir -p build
cd build

echo "==> Configuring build"
cmake .. \
    -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DREGRESSIONTEST_DOWNLOAD=OFF \
    -DCMAKE_BUILD_TYPE=Release

echo "==> Building"
make -j"$(nproc)"

echo "==> Installing"
make install

echo "==> Creating ~/bin/gmx launcher"
cat > "$BIN_DIR/gmx" <<EOF
#!/usr/bin/env bash
source "$PREFIX/bin/GMXRC"
exec "$PREFIX/bin/gmx" "\$@"
EOF
chmod +x "$BIN_DIR/gmx"

echo "Installation complete."
echo
echo "Then test with:"
echo "  gmx --version"
