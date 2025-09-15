#!/bin/bash
# Script for Julia Few Molecules spectrum calculation
# Author: Kacper Cybiński

# Local execution settings
WORKING_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/.."
cd "${WORKING_DIR}"

echo "========================================="
echo "Julia Few Molecules - Quick Local Test"
echo "Started at: " $(date)
echo "Working directory: " $(pwd)
echo "========================================="

# Set Julia multithreading
if [ -z "${JULIA_NUM_THREADS}" ]; then
    export JULIA_NUM_THREADS=11 # Number of threads fitted for M3 Pro CPU (11)
    echo "JULIA_NUM_THREADS set to $JULIA_NUM_THREADS"
else
    echo "JULIA_NUM_THREADS is already set to $JULIA_NUM_THREADS"
fi

# Global simulation parameters
MODE='2'        # Spectrum calculation
SUBTYPE='0'     # 0: vs Electric field, 1: vs Magnetic field, 2: vs A factor, 3: Geometry scaling, 4: vs Electric dipole, 5: vs Magnetic dipole
# Geometry: R_01 Φ_01 [R_02 Φ_02] [R_03 Φ_03] // R in nm, Φ in degrees
# GEOMETRY=''  # Single molecule at 0 nm, 0 degrees
GEOMETRY='200 90'  # Two molecules: one at 0 nm, 0 degrees, one at 500 nm, 90 degrees
# GEOMETRY='500 30 500 90'  # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, one at 30 degrees, one at 90 degrees
# EL_FIELD='0'  # Electric field in kV/m
# EL_FIELD='10'  # Electric field in kV/m
# EL_FIELD='11.40'  # Electric field in kV/m
# EL_FIELD='11.48'  # Electric field in kV/m
EL_FIELD='11.495'  # Electric field in kV/m
# EL_FIELD='50'  # Electric field in kV/m
PHI_EL='0'
THETA_EL='90'
# MG_FIELD='10'   # 10 Gauss magnetic field
MG_FIELD='50'   # 10 Gauss magnetic field
# MG_FIELD='0'   # 0 Gauss magnetic field for quick test
PHI_MG='0'
THETA_MG='90'
A_FACTOR='100'    # Spin-rotation coupling
# A_FACTOR='250'    # Spin-rotation coupling
D_EL='0.1'        # Electric dipole
D_MG='10'         # Magnetic dipole
M_TOT='-4.0'       # Total angular momentum projection
# M_TOT='1'       # Total angular momentum projection
START_IT='1'
END_IT='1001'      # Single iteration
# END_IT='501'      # Single iteration
# STEP='0.1'      # Small step for quick test
# STEP='0.005'      # Small step for quick test
# STEP='0.0002'      # Small step for quick test
# STEP='0.00004'      # Small step for quick test
# STEP='0.00003'      # Small step for quick test
STEP='0.00001'      # Small step for quick test
JMAX='2'       # Default from C++ constants.h
SPIN='5'      # Default from C++ constants.h
BROT='0.29'    # Not used in this test
DATETIME_NOW="$(date +"%d:%m:%Y_%H:%M:%S")"
ID="S$SPIN"
DTYPE="Float64"
PARALLEL_IT_ID=1001

# Check Julia installation
if ! command -v julia &> /dev/null; then
    echo "ERROR: Julia is not installed or not in PATH"
    echo "Please install Julia or add it to your PATH"
    exit 1
fi

# Check project structure
if [ ! -f "src/FewMolecules.jl" ]; then
    echo "ERROR: src/FewMolecules.jl not found. Make sure you're in the Julia_Few_Molecules directory."
    exit 1
fi

echo "Running quick test calculation..."
echo "Parameters: Mode=$MODE, Subtype=$SUBTYPE, B=$MG_FIELD G, geometry=$GEOMETRY"

mkdir -p logs

julia --project=. src/FewMolecules.jl \
    $MODE.$SUBTYPE \
    $GEOMETRY \
    --E $EL_FIELD \
    --phi_el $PHI_EL \
    --theta_el $THETA_EL \
    --B $MG_FIELD \
    --phi_mg $PHI_MG \
    --theta_mg $THETA_MG \
    --A $A_FACTOR \
    --d_el $D_EL \
    --d_mg $D_MG \
    --Mtot $M_TOT \
    --START_IT $START_IT \
    --END_IT $END_IT \
    --step $STEP \
    --JMAX $JMAX \
    --spin $SPIN \
    --basis-method "cpp-faithful" \
    --hamiltonian-method "cpp-faithful" \
    --hamiltonian-dtype $DTYPE \
    --looped-eigenproblem \
    --run-id-string $ID \
    --progressbar \
    --colored-print
job_status=$?
    # --IT_NO_OVERRIDE $PARALLEL_IT_ID \
# --run-id-string $DATETIME_NOW

echo "Test finished at: " $(date)
echo "Exit status: " $job_status
echo "========================================="

exit $job_status
