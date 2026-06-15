#!/bin/bash
set -e

TINKER_REPO="https://github.com/TinkerTools/tinker.git"
SRC_DIR="/tmp/tinker-build"
INSTALL_PREFIX="/opt/tinker"
PROFILE_FILE="$HOME/.tinker_env"

echo "==> Installing dependencies"
sudo apt-get update
sudo apt-get install -y \
  git \
  build-essential \
  gfortran \
  make \
  libfftw3-dev

echo "==> Cleaning old build dir"
rm -rf "$SRC_DIR"

echo "==> Cloning Tinker"
git clone "$TINKER_REPO" "$SRC_DIR"

echo "==> Preparing Makefile build"
cd "$SRC_DIR"

if [ ! -f "$SRC_DIR/make/Makefile" ]; then
  echo "ERROR: expected $SRC_DIR/make/Makefile but it was not found."
  exit 1
fi

if [ ! -d "$SRC_DIR/source" ]; then
  echo "ERROR: expected $SRC_DIR/source directory but it was not found."
  exit 1
fi

cp "$SRC_DIR/make/Makefile" "$SRC_DIR/source/Makefile"

echo "==> Patching Makefile paths"
sed -i "s|^TINKERDIR *=.*|TINKERDIR = $INSTALL_PREFIX|" "$SRC_DIR/source/Makefile"
sed -i "s|^LINKDIR *=.*|LINKDIR = $INSTALL_PREFIX/bin|" "$SRC_DIR/source/Makefile"

echo "==> Building Tinker"
cd "$SRC_DIR/source"
make -j"$(nproc)"

echo "==> Installing to $INSTALL_PREFIX"
sudo rm -rf "$INSTALL_PREFIX"
sudo mkdir -p "$INSTALL_PREFIX"
sudo cp -a "$SRC_DIR"/. "$INSTALL_PREFIX"/
sudo mkdir -p "$INSTALL_PREFIX/bin"

echo "==> Collecting executables"
find "$INSTALL_PREFIX/source" -maxdepth 1 -type f -name "*.x" -exec bash -c '
  for f in "$@"; do
    base="$(basename "$f" .x)"
    cp -f "$f" "'"$INSTALL_PREFIX"'/bin/$base"
    chmod +x "'"$INSTALL_PREFIX"'/bin/$base"
  done
' _ {} +

echo "==> Symlinking executables into /usr/local/bin"
for exe in "$INSTALL_PREFIX"/bin/*; do
  [ -f "$exe" ] || continue
  sudo ln -sf "$exe" "/usr/local/bin/$(basename "$exe")"
done

echo "==> Writing persistent environment file"
cat > "$PROFILE_FILE" <<EOF
export TINKER_HOME="$INSTALL_PREFIX"
export PATH="\$PATH:$INSTALL_PREFIX/bin"
EOF

if ! grep -q 'source ~/.tinker_env' "$HOME/.bashrc" 2>/dev/null; then
  echo '' >> "$HOME/.bashrc"
  echo '# Tinker environment' >> "$HOME/.bashrc"
  echo 'source ~/.tinker_env' >> "$HOME/.bashrc"
fi

echo "==> Loading environment for current shell"
# shellcheck disable=SC1090
source "$PROFILE_FILE"

echo "==> Done"
echo "TINKER_HOME=$TINKER_HOME"
echo "Try: analyze --help || ls $INSTALL_PREFIX/bin | head"
