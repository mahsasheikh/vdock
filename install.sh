#!/usr/bin/env bash
# Installs the "vina" environment for vdock. Linux, macOS (Intel + Apple Silicon), or WSL on Windows.
# Usage:  bash install.sh            (env name: vina)
#         ENV_NAME=vina127 bash install.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${ENV_NAME:-vina}"

case "$(uname -s)-$(uname -m)" in
  Darwin-arm64) ENV_FILE="$HERE/env/environment-osx-arm64.yml" ;;
  Linux-x86_64|Linux-aarch64|Darwin-x86_64) ENV_FILE="$HERE/env/environment.yml" ;;
  MINGW*|MSYS*|CYGWIN*) echo "Native Windows is not supported by conda-forge Vina. Use WSL2 (Ubuntu) and rerun there."; exit 1 ;;
  *) ENV_FILE="$HERE/env/environment.yml" ;;
esac

# Pick a package manager: micromamba > mamba > conda. Install micromamba if none.
if command -v micromamba >/dev/null 2>&1; then PM=micromamba
elif command -v mamba >/dev/null 2>&1; then PM=mamba
elif command -v conda >/dev/null 2>&1; then PM=conda
else
  echo ">> No conda-type manager found; installing micromamba (official installer)."
  "${SHELL:-bash}" <(curl -L micro.mamba.pm/install.sh) </dev/null
  export PATH="$HOME/.local/bin:$PATH"
  PM=micromamba
fi
echo ">> Using $PM with $ENV_FILE -> env '$ENV_NAME'"

if $PM env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo ">> Env '$ENV_NAME' exists; updating."
  $PM env update -y -n "$ENV_NAME" -f "$ENV_FILE" 2>/dev/null || $PM install -y -n "$ENV_NAME" -f "$ENV_FILE"
else
  $PM env create -y -n "$ENV_NAME" -f "$ENV_FILE" 2>/dev/null || $PM create -y -n "$ENV_NAME" -f "$ENV_FILE"
fi

echo ">> Verifying installation"
$PM run -n "$ENV_NAME" bash -c '
  set -e
  vina --version
  for t in mk_prepare_receptor.py mk_prepare_ligand.py mk_export.py scrub.py autogrid4; do
    command -v $t >/dev/null || { echo "MISSING: $t"; exit 1; }; done
  python -c "import vina, meeko, rdkit, molscrub, pdbfixer, gemmi, prolif; print(\"vina\", vina.__version__, \"| meeko\", meeko.__version__, \"| rdkit\", rdkit.__version__)"
'
$PM run -n "$ENV_NAME" python -m pip freeze > "$HERE/env/installed-$(uname -s)-$(uname -m).lock.txt" 2>/dev/null || true
echo
echo "Installed. Next:"
echo "  $PM activate $ENV_NAME"
echo "  bash $HERE/examples/run_1iep_redock.sh      # ~3-10 min validation"
