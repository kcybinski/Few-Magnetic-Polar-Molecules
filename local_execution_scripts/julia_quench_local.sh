#!/bin/bash
# Local Julia Few Molecules quench dynamics script
# Author: Kacper Cybiński

# Local execution settings
WORKING_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/.."
cd "${WORKING_DIR}"

# Local parameters for quench dynamics
MODE=1  # 1 - quench dynamics
SUBTYPE=0  # Quench type: 0:ms and m, 1: Ms and M
# GEOMETRY=''  # geometry: R_01 Φ_01 [R_02 Φ_02] [R_03 Φ_03]
GEOMETRY='100 0'  # geometry: R_01 Φ_01 [R_02 Φ_02] [R_03 Φ_03]
# GEOMETRY='500 30 500 90'  # geometry: R_01 Φ_01 [R_02 Φ_02] [R_03 Φ_03]
ELF_ANGLES='0 90'  # phi_el theta_el
MF_ANGLES='0 90'  # phi_mg theta_mg
RRR=500

JMAX=2  # Maximum total angular momentum (default from C++ constants.h)
SPIN=5  # Spin value (default from C++ constants.h)
# MTOT=-0.5  # Total angular momentum projection (default from C++ constants.h)

BROT=0.29  # Spin-rotation coupling constant in GHz (default from C++ constants.h, not used in scripts)
HARMONIC_FREQ=100.0  # Trap harmonic frequency in kHz (default, can be overridden)
HARMONIC_FREQ_UNIT="kHz"  # Unit for harmonic frequency (default, can be overridden)

D_EL=0.1  # Electric dipole moment in Debye
D_MG=10  # Magnetic dipole moment in μB
A_FACTOR=100  # Spin-rotation coupling in kHz
MG_FIELD=50  # Magnetic field in Gauss
# EL_FIELD=11.49125  # Electric field strength in kV/cm
# EL_FIELD=11.49125  # Electric field strength in kV/cm
EL_FIELD=11.4848  # Electric field strength in kV/cm

quenched_B=50  # Quenched magnetic field in Gauss (default, can be overridden)
# quenched_E=11.49251  # Quenched electric field in kV/cm (default, can be overridden)
# quenched_E=11.49293  # Quenched electric field in kV/cm (default, can be overridden)
quenched_E=11.49914  # Quenched electric field in kV/cm (default, can be overridden)
quenched_R='200 90'  # Quenched geometry radius in nm (default, can be overridden)

# 2 mol
# product_state="1,-1,2.5,-2.5,1,1,2.5,2.5"  # Product state for quench dynamics --> State (1)
# MTOT=0.0
# product_state="1,1,2.5,2.5,1,-1,2.5,0.5"  # Product state for quench dynamics --> State (2)
# MTOT=3.0
product_state="1,-1,2.5,-2.5,1,-1,2.5,0.5"  # Product state for quench dynamics --> State (3)
MTOT=-4.0

# 1 mol
# product_state="1,-1,2.5,-2.5"  # Product state for quench dynamics --> Mol (1)
# MTOT=-3.5
# product_state="1,1,2.5,2.5"  # Product state for quench dynamics --> Mol (2)
# MTOT=3.5
# product_state="1,-1,2.5,0.5"  # Product state for quench dynamics --> Mol (3)
# MTOT=-0.5


# eigenstate_no=1365  # Initial eigenstate number 
eigenstate_no=76 # Initial eigenstate number (1-based indexing)
basis_state_no=59 # Initial basis state number (1-based indexing)

# Time parameters for quench (in harmonic units)
TMAX=5000.0    # Maximum time in harmonic units
DT=0.25        # Time step in harmonic units
T_UNIT="harmonic"  # Time unit for quench dynamics

# # (Alternative) Time parameters in SI units
# TMAX=25.0  # Maximum time
# T_UNIT="μs" # "Input time unit: 'harmonic' (dimensionless 2π/ω units) or SI units ('s', 'ms', 'μs', 'ns', 'ps', 'fs')"
# DT=0.01  # Time step in SI units

# Set Julia multithreading
if [ -z "${JULIA_NUM_THREADS}" ]; then
    export JULIA_NUM_THREADS=11  # Conservative for test
fi

echo "========================================="
echo "Julia Few Molecules - Local Quench Dynamics"
echo "Started at: " $(date)
echo "Working directory: " $(pwd)
echo "WARNING: Quench dynamics can be computationally intensive"
echo "========================================="

echo "Starting quench dynamics calculations..."

# Generate unique job name
JOBNAME="QT${MODE}.${SUBTYPE}DEL${D_EL}D${A_FACTOR}M${M_TOT}DMG${D_MG}B${MG_FIELD}E${EL_FIELD}R${RRR}"

echo "[$job_count/$total_jobs] Running quench dynamics job: $JOBNAME"

echo "Starting quench dynamics with: D_EL=$D_EL, D_MG=$D_MG, A=$A_FACTOR, B=$MG_FIELD, E=$EL_FIELD..."
if [ "$T_UNIT" == "harmonic" ]; then
    echo "Using harmonic units for time: TMAX=$TMAX × (2π/ω), DT=$DT × (2π/ω)"
else
    echo "Using SI units for time: TMAX=$TMAX $T_UNIT, DT=$DT $T_UNIT"
fi

# Run Julia quench dynamics with CLI arguments only
julia --project=. src/FewMolecules.jl \
    $MODE.$SUBTYPE \
    $GEOMETRY \
    --E $EL_FIELD \
    --phi_el $(echo ${ELF_ANGLES} | cut -d' ' -f1) \
    --theta_el $(echo ${ELF_ANGLES} | cut -d' ' -f2) \
    --B $MG_FIELD \
    --phi_mg $(echo ${MF_ANGLES} | cut -d' ' -f1) \
    --theta_mg $(echo ${MF_ANGLES} | cut -d' ' -f2) \
    --A $A_FACTOR \
    --d_el $D_EL \
    --d_mg $D_MG \
    --Mtot $MTOT \
    --JMAX $JMAX \
    --spin $SPIN \
    --Brot $BROT \
    --basis-method cpp-faithful \
    --hamiltonian-method cpp-faithful \
    --eigensolver auto \
    --hamiltonian-dtype "ComplexF64" \
    --harmonic-frequency $HARMONIC_FREQ \
    --harmonic-frequency-unit $HARMONIC_FREQ_UNIT \
    --harmonic-frequency-source custom \
    --tmax $TMAX \
    --dt $DT \
    --quenched_E $quenched_E \
    --quenched_R $quenched_R \
    --productstate "$product_state" \
    --quench-observables M Ms ms1 m1 ms2 m2\
    --time-unit "harmonic" \
    --progressbar \
    # --initial-eigenstate $eigenstate_no \
    # --initial-basis-state $basis_state_no \
    # --quench-observables M Ms ms1 m1\
    # --quench-observables M Ms ms1 m1 ms2 m2 ms3 m3\
    # --quenched_R $quenched_R \
    # --skip-terms "electric_dipolar, magnetic_dipolar" \
    # Uncomment the next line to enable profiling
    # --profile

if [ $? -eq 0 ]; then
    echo "✓ Quench dynamics job $JOBNAME completed successfully"
else
    echo "✗ Quench dynamics job $JOBNAME failed"
fi

echo "---"

echo "========================================="
echo "All quench dynamics jobs completed!"
echo "Finished at: " $(date)
echo "Results stored in data/quench_results/"
echo "Check individual output files for detailed results"
echo "========================================="
echo "Results stored in data/quench_results/"
echo "Check individual output files for detailed results"
echo "========================================="
echo "========================================="
