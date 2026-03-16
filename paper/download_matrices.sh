#!/bin/bash
# Download test matrices from SuiteSparse Matrix Collection
# Usage: ./download_matrices.sh [output_dir]

set -e

OUTDIR="${1:-matrices}"
mkdir -p "$OUTDIR"

BASE="https://sparse.tamu.edu/MM"

# All test matrices used in the propack++ paper benchmarks.
# Format: "Name Group/Name"
MATRICES=(
  # Original PROPACK test matrices (small)
  "well1850 HB/well1850"
  "illc1850 HB/illc1850"
  "tols4000 Bai/tols4000"
  "mhd4800a Bai/mhd4800a"
  # Medium matrices
  "rdb5000 Bai/rdb5000"
  "bcsstk38 Boeing/bcsstk38"
  "dw8192 Bai/dw8192"
  "ex19 FIDAP/ex19"
  "af23560 Bai/af23560"
  "email-Enron SNAP/email-Enron"
  "cant Williams/cant"
  "rail_79841 Oberwolfach/rail_79841"
  # Large square unsymmetric
  "torso1 Norris/torso1"
  "stomach Norris/stomach"
  # Large square symmetric (structural/FEM)
  "inline_1 GHS_psdef/inline_1"
  "Fault_639 Janna/Fault_639"
  "Emilia_923 Janna/Emilia_923"
  "bone010 Oberwolfach/bone010"
  # Large rectangular (standard SVD benchmarks)
  "LargeRegFile Stevenson/LargeRegFile"
  "Rucci1 Rucci/Rucci1"
)

download_matrix() {
  local name="$1"
  local group_name="$2"

  if [ -f "$OUTDIR/${name}.mtx" ]; then
    echo "  $name.mtx already exists, skipping."
    return 0
  fi

  local url="$BASE/${group_name}.tar.gz"
  local tarball="/tmp/${name}.tar.gz"

  echo "  Downloading $name from $url ..."
  if ! curl -sS -L -o "$tarball" "$url"; then
    echo "  ERROR: Failed to download $name"
    return 1
  fi

  # Check for HTML error page
  if file "$tarball" | grep -q "HTML"; then
    echo "  ERROR: Got HTML instead of tarball for $name (404?)"
    rm -f "$tarball"
    return 1
  fi

  echo "  Extracting..."
  tar -xzf "$tarball" -C "$OUTDIR" 2>/dev/null || true

  # Find the .mtx file
  local mtx_file
  mtx_file=$(find "$OUTDIR" -name "${name}.mtx" -type f 2>/dev/null | head -1)
  if [ -z "$mtx_file" ]; then
    mtx_file=$(find "$OUTDIR" -iname "*.mtx" -path "*${name}*" -type f 2>/dev/null | head -1)
  fi

  if [ -n "$mtx_file" ] && [ "$mtx_file" != "$OUTDIR/${name}.mtx" ]; then
    mv "$mtx_file" "$OUTDIR/${name}.mtx"
  fi

  # Clean up
  rm -f "$tarball"
  find "$OUTDIR" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} + 2>/dev/null || true

  if [ -f "$OUTDIR/${name}.mtx" ]; then
    local size
    size=$(du -h "$OUTDIR/${name}.mtx" | cut -f1)
    echo "  OK: ${name}.mtx ($size)"
  else
    echo "  WARNING: ${name}.mtx not found after extraction"
  fi
}

echo "Downloading matrices to $OUTDIR/"
echo "============================================"

for entry in "${MATRICES[@]}"; do
  read -r name group <<< "$entry"
  download_matrix "$name" "$group"
done

echo ""
echo "Done. Contents of $OUTDIR/:"
ls -lhS "$OUTDIR"/*.mtx 2>/dev/null || echo "  (no .mtx files found)"
