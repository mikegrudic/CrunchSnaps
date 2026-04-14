#!/bin/bash
# Regenerate all example images for the documentation.
# Run from the CrunchSnaps repo root with a snapshot_2000.hdf5 present.
set -e

SNAP=snapshot_2000.hdf5
RES=1024
OUT=docs/_static/examples
SINKVIS2="SinkVis2"

if [ ! -f "$SNAP" ]; then
    echo "Error: $SNAP not found in $(pwd)" >&2
    exit 1
fi

mkdir -p "$OUT"
export PYTHONHASHSEED=0

run() {
    local name="$1"; shift
    rm -rf "$OUT/.maps"
    echo "Generating $name..."
    "$SINKVIS2" "$@" 2>&1 | grep -v "^$"
    # Rename the output to the desired name
    local src=$(ls -t "$OUT"/*.png 2>/dev/null | head -1)
    if [ -n "$src" ] && [ "$src" != "$OUT/$name" ]; then
        mv "$src" "$OUT/$name"
    fi
}

# Default SigmaGas (viridis)
run sigmagas_default.png "$SNAP" --res=$RES --outputfolder="$OUT"

# Inferno colormap
run sigmagas_inferno.png "$SNAP" --res=$RES --cmap=inferno --outputfolder="$OUT"

# Slice(Temperature)
run slice_temperature.png "$SNAP" --tasks='Slice(Temperature)' --res=$RES --outputfolder="$OUT"

# Slice(PlasmaBeta) with RdBu_r
# NOTE: --unit_B=1e4 because this snapshot stores B in Tesla, not Gauss
run slice_plasmabeta.png "$SNAP" --tasks='Slice(PlasmaBeta)' --res=$RES --cmap=RdBu_r --unit_B=1e4 --outputfolder="$OUT"

# ProjectedAverage(Temperature)
run projectedavg_temperature.png "$SNAP" --tasks='ProjectedAverage(Temperature)' --res=$RES --outputfolder="$OUT"

# SurfaceDensity(Masses*InternalEnergy) with magma
run surfacedensity_thermal.png "$SNAP" --tasks='SurfaceDensity(Masses*InternalEnergy)' --res=$RES --cmap=magma --outputfolder="$OUT"

# View along x
run sigmagas_xdir.png "$SNAP" --res=$RES --dir=x --outputfolder="$OUT"

# Pan and tilt
run sigmagas_pantilt.png "$SNAP" --res=$RES --pan=45 --tilt=15 --outputfolder="$OUT"

# Center on densest, zoomed in
run sigmagas_densest.png "$SNAP" --res=$RES --center=densest --rmax=5 --outputfolder="$OUT"

# Slice(MachNumber) with inferno
# MachNumber doesn't depend on B units, but keep consistent
run slice_machnumber.png "$SNAP" --tasks='Slice(MachNumber)' --res=$RES --cmap=inferno --outputfolder="$OUT"

# Sigma1D (velocity dispersion)
run sigma1d_velocities.png "$SNAP" --tasks='Sigma1D(Velocities)' --res=$RES --cmap=magma --outputfolder="$OUT"

# Matplotlib backend
run sigmagas_matplotlib.png "$SNAP" --res=$RES --backend=matplotlib --outputfolder="$OUT"

# Slice(Entropy) with plasma
run slice_entropy.png "$SNAP" --tasks='Slice(Entropy)' --res=$RES --cmap=plasma --outputfolder="$OUT"

echo ""
echo "Generated $(ls "$OUT"/*.png | wc -l) example images in $OUT/"
ls "$OUT"/*.png
