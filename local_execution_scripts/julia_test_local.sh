#!/bin/bash
# Quick test script for Julia Few Molecules - local execution
# Simple single calculation to verify everything works
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
    export JULIA_NUM_THREADS=11  # Conservative for test
fi

# Test parameters - simple case
MODE='0'        # Spectrum calculation
SUBTYPE='0'     # vself
# GEOMETRY=''  # Simple geometry
GEOMETRY='100 90'  # Simple geometry
# GEOMETRY='500 30 500 90'  # Two molecules at 500 nm, one at 30 degrees, one at 90 degrees
EL_FIELD='0'  # Electric field in V/m
PHI_EL='0'
THETA_EL='0'
MG_FIELD='10'   # 10 Gauss magnetic field
PHI_MG='0'
THETA_MG='0'
A_FACTOR='100'    # Spin-rotation coupling
D_EL='0.1'        # Electric dipole
D_MG='10'        # Magnetic dipole
M_TOT='-1.0'       # Total angular momentum projection
START_IT='1'
END_IT='1'      # Single iteration
STEP='0.01'      # Small step for quick test
JMAX='2'       # Default from C++ constants.h
SPIN='5'      # Default from C++ constants.h
BROT='0.29'    # Not used in this test

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

# Create directories
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
    --hamiltonian-dtype "Float64" \
    --looped-eigenproblem \
    --progressbar \
    --colored-print \
    --profile
job_status=$?

echo "Test finished at: " $(date)
echo "Exit status: " $job_status
echo "========================================="

exit $job_status
