#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
MINIPROT_VERSION="v0.18"

python3 -m venv "$ROOT/.venv"
source "$ROOT/.venv/bin/activate"

python -m pip install --upgrade pip
python -m pip install lifton

mkdir -p "$ROOT/.tools" "$ROOT/.build"

if [ ! -x "$ROOT/.tools/miniprot" ]; then
    rm -rf "$ROOT/.build/miniprot"

    git clone \
        --depth 1 \
        --branch "$MINIPROT_VERSION" \
        https://github.com/lh3/miniprot.git \
        "$ROOT/.build/miniprot"

    make -C "$ROOT/.build/miniprot"
    cp "$ROOT/.build/miniprot/miniprot" "$ROOT/.tools/miniprot"
fi

echo
echo "Setup complete."
echo "Activate with:"
echo "  source $ROOT/.venv/bin/activate"
echo "  export PATH=\"$ROOT/.tools:\$PATH\""
