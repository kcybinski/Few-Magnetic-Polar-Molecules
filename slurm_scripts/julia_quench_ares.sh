#!/bin/bash
# SLURM script for Julia Few Molecules quench dynamics on PLGrid Ares
# Generates parameter sweep jobs using JuliaFewMolecules template
# Author: Kacper Cybiński

# Global simulation parameters for easy editing
MODE='1'  # 1 - quench dynamics
SUBTYPE='0'  # 0: calculate , 1: consolidate distributed results, 2: plot result time series and Fourier transforms
# Geometry: R_01 Φ_01 [R_02 Φ_02] [R_03 Φ_03] // R in nm, Φ in degrees
# GEOMETRY=''  # Single molecule at 0 nm, 0 degrees
GEOMETRY='500 90'  # Two molecules: one at 0 nm, 0 degrees, one at 500 nm, 90 degrees
# GEOMETRY='500 30 500 90'  # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, one at 30 degrees, one at 90 degrees
N_MOL=2

D_EL='0.1'        # Electric dipole
D_MG='10'         # Magnetic dipole
M_TOT='-1'       # Total angular momentum projection
JMAX='2'  # Maximum total angular momentum
SPIN='5'  # Spin value (2S)
BROT='0.29'  # Spin-rotation coupling constant in GHz

# # External fields
EL_FIELD='11.46'  # Electric field in kV/m
PHI_EL='0' # Electric field angles
THETA_EL='0' # Electric field angles

MG_FIELD='10'   # 10 Gauss magnetic field
PHI_MG='0' # Magnetic field angles
THETA_MG='0' # Magnetic field angles

A_FACTOR='100'    # Spin-rotation coupling

HARMONIC_FREQ='100.0'  # Trap harmonic frequency in kHz (default, can be overridden)
HARMONIC_UNIT='kHz'  # Unit for harmonic frequency

# Quench parameters
# quenched_B_val=10  # Quenched magnetic field in Gauss (default, can be overridden)
quenched_E_val=11.50  # Quenched electric field in kV/cm (default, can be overridden)
# quenched_R_val=$(echo $GEOMETRY | awk '{print $1}')  # Quenched geometry radius in nm (default, can be overridden)

# State specification
# product_state="1,0,6.5,-6.5,1,-1,6.5,6.5"  # Product state for quench dynamics
# job_state_id="PS${product_state}"

# eigenstate_no=1365  # Initial eigenstate number 
eigenstate_no=189  # Initial eigenstate number 
job_state_id="ES${eigenstate_no}"


# Time parameters for quench (in harmonic units)
TMAXval=80.0    # Maximum time in harmonic units
DT_VAL=2.0        # Time step in harmonic units
T_UNIT="harmonic"  # Time unit for quench dynamics
DTYPE="ComplexF64" # Format of the Hamiltonian matrix

OBSERVABLES_LIST="M Ms ms1 m1 ms2 m2 ms3 m3"

# SLURM job parameters
NRAM='500'  # [MB]
TIME='00:01:00'  # format HH:MM:SS
QUEUE='plgrid'
GRANT='plgquantmol8-cpu'
NCORES=3 # Number of cores/threads to use
DTNOW="$(date +"%d-%m-%Y_%H-%M-%S")"
IDSTR="S$SPIN"

echo "========================================="
echo "Julia Few Molecules - SLURM Quench Dynamics Setup"
echo "Started at: " $(date)
echo "Setting up quench dynamics jobs for PLGrid Ares"
echo "WARNING: Quench dynamics can be computationally intensive"
echo "========================================="

# Create output directories if they don't exist
mkdir -p ../logs
mkdir -p ../logs/outputs
mkdir -p ../logs/outputs/DTNOW
mkdir -p ../logs/errors
mkdir -p ../logs/errors/DTNOW

# Generate unique job name
JOBNAME="QD${MODE}.${SUBTYPE}_Jmax${JMAX}_Brot${BROT}_s${SPIN}_E${EL_FIELD}_B${MG_FIELD}_A${A_FACTOR}_del${D_EL}_dmg${D_MG}_Mtot${M_TOT}_NMOL${N_MOL}_${job_state_id}"

# Create job script from template
sed -e "s/JOBNAME/${JOBNAME}/g" \
    -e "s/TIME/${TIME}/g" \
    -e "s/RAM/${NRAM}/g" \
    -e "s/GRANTID/${GRANT}/g" \
    -e "s/QUEUEID/${QUEUE}/g" \
    -e "s/NCORES/${NCORES}/g" \
    -e "s/MODE.SUBTYPE/${MODE}.${SUBTYPE}/g" \
    -e "s/GEOMETRY/${GEOMETRY}/g" \
    -e "s/EL_FIELD/${EL_FIELD}/g" \
    -e "s/PHI_EL/${PHI_EL}/g" \
    -e "s/THETA_EL/${THETA_EL}/g" \
    -e "s/MG_FIELD/${MG_FIELD}/g" \
    -e "s/PHI_MG/${PHI_MG}/g" \
    -e "s/THETA_MG/${THETA_MG}/g" \
    -e "s/A_FACTOR/${A_FACTOR}/g" \
    -e "s/D_EL/${D_EL}/g" \
    -e "s/D_MG/${D_MG}/g" \
    -e "s/MTOT/${M_TOT}/g" \
    -e "s/JMAXval/${JMAX}/g" \
    -e "s/SPINval/${SPIN}/g" \
    -e "s/DTYPE/${DTYPE}/g" \
    -e "s/HARMONIC_FREQ/${HARMONIC_FREQ}/g" \
    -e "s/HARMONIC_UNIT/${HARMONIC_UNIT}/g" \
    -e "s/T_UNIT/${T_UNIT}/g" \
    -e "s/TMAXval/${TMAXval}/g" \
    -e "s/DT_VAL/${DT_VAL}/g" \
    -e "s/quenched_B_val/${quenched_B_val}/g" \
    -e "s/quenched_E_val/${quenched_E_val}/g" \
    -e "s/quenched_R_val/${quenched_R_val}/g" \
    -e "s/OBSERVABLES_LIST/${OBSERVABLES_LIST}/g" \
    -e "s/EIGENSTATE_NO/${eigenstate_no}/g" \
    -e "s/PRODUCT_STATE/${product_state}/g" \
    -e "s/IDSTR/${IDSTR}/g" \
    -e "s/DTNOW/${DTNOW}/g" \
    JuliaFewMoleculesQuench > ${JOBNAME}.slurm

# Submit job
sbatch ${JOBNAME}.slurm

if [ $? -eq 0 ]; then
    echo "✓ Job $JOBNAME submitted successfully"
else
    echo "✗ Failed to submit job $JOBNAME"
fi

# Clean up job script
rm -f ${JOBNAME}.slurm

# Brief pause to avoid overwhelming the scheduler
sleep 0.1

echo "========================================="
echo "All quench dynamics jobs submitted!"
echo "Finished at: " $(date)
echo "Total jobs submitted: $job_count"
echo "Monitor jobs with: squeue -u \$USER"
echo "Results will be in ${SCRATCH}/Julia_Few_Molecules//data/spectrum_results/ with organized hierarchy"
echo "========================================="
