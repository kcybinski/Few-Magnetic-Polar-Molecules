#!/bin/bash
# SLURM script for Julia Few Molecules spectrum calculation on PLGrid Ares
# Generates parameter sweep jobs using JuliaFewMolecules template
# Author: Kacper Cybiński

# Global simulation parameters for easy editing
MODE='2'        # 0 - spectrum calculation (energies), 2 - spectrum and eigenvalues calculation
SUBTYPE='0'     # 0:vself, 1:vsmf, 2:vsA, 3:vs_geometry, 4:vsd_el, 5:vsd_mg

# # Molecule parameters

# Geometry: R_01 Φ_01 [R_02 Φ_02] [R_03 Φ_03] // R in nm, Φ in degrees
# GEOMETRY=''  # Single molecule at 0 nm, 0 degrees
# GEOMETRY='500 90'  # Two molecules: one at 0 nm, 0 degrees, one at 500 nm, 90 degrees
# GEOMETRY='500 30 500 90'  # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, 30 degrees, one at 500 nm, 90 degrees (equilateral triangle geometry)
GEOMETRY='100 0 100 180'  # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, 0 degrees, one at 1000 nm, 180 degrees (planar vertical geometry)
N_MOL=3

D_EL='0.1'        # Electric dipole
D_MG='10'         # Magnetic dipole
M_TOT='-0.5'       # Total angular momentum projection
JMAX='2'  # Maximum total angular momentum
SPIN='5'  # Spin value (2S)
BROT='0.29'  # Spin-rotation coupling constant in GHz

# # External fields
EL_FIELD='0'  #  Electric field in kV/m
# EL_FIELD='11.38'  # Electric field in kV/m
PHI_EL='0' # Electric field angles
THETA_EL='90' # Electric field angles

MG_FIELD='50'   # 50 Gauss magnetic field
PHI_MG='0' # Magnetic field angles
THETA_MG='90' # Magnetic field angles

A_FACTOR='100'    # Spin-rotation coupling

# # Simulation  parameters
START_IT=1 # >= 1, determines from where looped calculation runs
# END_IT=2 # Maximal iteration index
END_IT=501 # Maximal iteration index
STEP=0.03 # step size parameter for looped calculations
DTYPE="Float64" # Format of the Hamiltonian matrix
NUM_PTS_IN_LOOP_PER_THREAD=1 # Number of points in the loop per thread. This allows for calculating more slices in the diagram than the 1000 jobs cap on Ares


# SLURM job parameters
# NRAM='2500'  # [MB]
# TIME='00:02:00'  # format HH:MM:SS
# NCORES=3 # Number of cores/threads to use

# Spectrum calculation for 3 molecules, s=5, mode = 2
NRAM='14000'  # [MB]
TIME='00:25:00'  # format HH:MM:SS
NCORES=3 # Number of cores/threads to use

# Spectrum calculation for 3 molecules, s=8, mode = 2
# NRAM='50000'  # [MB]
# TIME='06:00:00'  # format HH:MM:SS
# NCORES=3 # Number of cores/threads to use

# # Code testing parameters
# NRAM='1000'  # [MB]
# TIME='00:01:00'  # format HH:MM:SS
# NCORES=1 # Number of cores/threads to use

QUEUE='YOUR_QUEUE'  # e.g. 'plgrid'
GRANT='YOUR_GRANT_NAME'
DTNOW="$(date +"%d-%m-%Y_%H-%M-%S")"
IDSTR="S$SPIN"

total_jobs=$(END_IT-START_IT)

echo "========================================="
echo "Julia Few Molecules - SLURM Spectrum Calculation Setup"
echo "Started at: " $(date)
echo "Setting up parameter sweep jobs for PLGrid Ares"
echo "========================================="

echo "Total spectrum calculation jobs to submit: $total_jobs"
echo "Generating and submitting SLURM jobs..."

# Create output directories if they don't exist
mkdir -p ../logs
mkdir -p ../logs/outputs
mkdir -p ../logs/outputs/DTNOW
mkdir -p ../logs/errors
mkdir -p ../logs/errors/DTNOW


# The total number of jobs is the number of iterations by iterations per job
if [ "$NUM_PTS_IN_LOOP_PER_THREAD" -eq 1 ]; then
    total_jobs=$(printf "%.0f" $(echo "$END_IT - $START_IT" | bc))
else
    total_jobs=$(printf "%.0f" $(echo "(($END_IT - $START_IT) / $NUM_PTS_IN_LOOP_PER_THREAD) - 1" | bc))
fi

job_count=0
for PARALLEL_IT_ID in $(seq 0 $total_jobs); do
    # Generate unique job name
    # JOBNAME="SP${MODE}.${SUBTYPE}DEL${D_EL}D${A_FACTOR}M${M_TOT}DMG${D_MG}B${MG_FIELD}S${START_IT}E${EL_FIELD}"


    JOBNAME="${job_count}_SP${MODE}.${SUBTYPE}_Jmax${JMAX}_Brot${BROT}_s${SPIN}_E${EL_FIELD}_B${MG_FIELD}_A${A_FACTOR}_del${D_EL}_dmg${D_MG}_Mtot${M_TOT}_NMOL${N_MOL}_GEOMRq_$(echo "$GEOMETRY" | tr ' ' '_')"

    echo "[$job_count/$total_jobs] Creating SLURM job: $JOBNAME"
    job_count=$((job_count + 1))

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
        -e "s/START_IT_VALUE/${START_IT}/g" \
        -e "s/END_IT_VALUE/${END_IT}/g" \
        -e "s/STEP/${STEP}/g" \
        -e "s/JMAXval/${JMAX}/g" \
        -e "s/SPINval/${SPIN}/g" \
        -e "s/DTYPE/${DTYPE}/g" \
        -e "s/IDSTR/${IDSTR}/g" \
        -e "s/DTNOW/${DTNOW}/g" \
        -e "s/ITNUMVAL/${PARALLEL_IT_ID}/g" \
       -e "s/ITSNUMOVERRIDE/${NUM_PTS_IN_LOOP_PER_THREAD}/g" \
        JuliaFewMoleculesDistributed > ${JOBNAME}.slurm
    
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
    sleep 0.05
done

echo "========================================="
echo "All spectrum calculation jobs submitted!"
echo "Finished at: " $(date)
echo "Total jobs submitted: $job_count"
echo "Monitor jobs with: hpc-jobs"
echo "Results will be in ${SCRATCH}/Julia_Few_Molecules/data/spectrum_resuls/ with organized hierarchy"
echo "========================================="
