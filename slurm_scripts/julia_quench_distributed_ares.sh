#!/bin/bash
# SLURM script for Julia Few Molecules quench dynamics on PLGrid Ares
# Generates parameter sweep jobs using JuliaFewMolecules template
# Author: Kacper Cybiński

# Global simulation parameters for easy editing
MODE='1'  # 1 - quench dynamics
SUBTYPE='0' # 0: calculate , 1: consolidate distributed results, 2: plot result time series and Fourier transforms
# Geometry: R_01 Φ_01 [R_02 Φ_02] [R_03 Φ_03] // R in nm, Φ in degrees
# GEOMETRY=''  # Single molecule at 0 nm, 0 degrees
# GEOMETRY='500 90'  # Two molecules: one at 0 nm, 0 degrees, one at 500 nm, 90 degrees
GEOMETRY='500 30 500 90'  # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, one at 30 degrees, one at 90 degrees
N_MOL=3

D_EL='0.1'        # Electric dipole
D_MG='10'         # Magnetic dipole
M_TOT='-0.5'       # Total angular momentum projection
JMAX='2'  # Maximum total angular momentum
SPIN='5'  # Spin value (2S)
BROT='0.29'  # Spin-rotation coupling constant in GHz

# # External fields
EL_FIELD='11.47292'  # Electric field in kV/m
PHI_EL='0' # Electric field angles
THETA_EL='90' # Electric field angles

MG_FIELD='50'   # 10 Gauss magnetic field
PHI_MG='0' # Magnetic field angles
THETA_MG='90' # Magnetic field angles

A_FACTOR='100'    # Spin-rotation coupling

HARMONIC_FREQ='100.0'  # Trap harmonic frequency in kHz (default, can be overridden)
HARMONIC_UNIT='kHz'  # Unit for harmonic frequency

# Quench parameters
# quenched_B_val=10  # Quenched magnetic field in Gauss (default, can be overridden)
quenched_E_val=11.4927 # Quenched electric field in kV/cm (default, can be overridden)
# quenched_E_val=11.49914 # Quenched electric field in kV/cm (default, can be overridden)
quenched_R_val='100 0 100 180'  # Quenched geometry radius in nm (default, can be overridden)

# State specification
# product_state="1,-1,2.5, -2.5, 1, 1, 2.5, 2.5, 1, -1, 2.5, 0.5"  # Product state for quench dynamics
# # job_state_id="PS${product_state}"
# job_state_id="PS"

# eigenstate_no=1883  # Initial eigenstate number for 3 molecules
# eigenstate_no=189  # Initial eigenstate number for 2 molecules 
# job_state_id="ES${eigenstate_no}"

basis_state_no=3035  # Initial basis state number 
job_state_id="BS${basis_state_no}"


# Time parameters for quench (in harmonic units)
TMAXval=250.0    # Maximum time in harmonic units
DT_VAL=0.25        # Time step in harmonic units
NUM_PTS_IN_LOOP_PER_THREAD=10 # Number of points in the loop per thread. This allows for calculating more slices in the diagram than the 1000 jobs cap on Ares

# Values for testing code validity
# TMAXval=2.0    # Maximum time in harmonic units
# DT_VAL=0.5        # Time step in harmonic units
# NUM_PTS_IN_LOOP_PER_THREAD=2 # Number of points in the loop per thread. This allows for calculating more slices in the diagram than the 1000 jobs cap on Ares

T_UNIT="harmonic"  # Time unit for quench dynamics
DTYPE="ComplexF64" # Format of the Hamiltonian matrix

OBSERVABLES_LIST="M Ms ms1 m1 ms2 m2 ms3 m3" # Three molecules 
# OBSERVABLES_LIST="M Ms ms1 m1 ms2 m2" # Two molecules

# SLURM job parameters
# NRAM='2500'  # [MB]
# TIME='00:03:00'  # format HH:MM:SS

# Values for S=5 quench
NRAM='60000'  # [MB]
TIME='04:00:00'  # format HH:MM:SS
NCORES=3 # Number of cores/threads to use

# Code test values
# NRAM='2000'  # [MB]
# TIME='00:05:00'  # format HH:MM:SS
# NCORES=1 # Number of cores/threads to use

DTNOW="$(date +"%d-%m-%Y_%H-%M-%S")"
IDSTR="S$SPIN"

# Parameters for SUBTYPE=1, consolidation of results
# NRAM_CONSOLIDATE='2000'
# TIME="00:03:00"
# NCORES=1 # Single core for consolidation

# Constant grant parameters
QUEUE='YOUR_QUEUE'  # e.g. 'plgrid'
GRANT='YOUR_GRANT_NAME'

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

# The total number of jobs is the number of time steps, from 0 to TMAXval
if [ "$NUM_PTS_IN_LOOP_PER_THREAD" -eq 1 ]; then
    total_jobs=$(echo "($TMAXval / $DT_VAL)" | bc)
else
    total_jobs=$(echo "($TMAXval / ($DT_VAL * $NUM_PTS_IN_LOOP_PER_THREAD)) - 1" | bc)
fi


job_count=0

if [ "$SUBTYPE" -eq 0 ]; then
    for PARALLEL_IT_ID in $(seq 0 $total_jobs); do

        IT_ID_NOW_INT=$(printf "%.0f" $(echo "$PARALLEL_IT_ID * $NUM_PTS_IN_LOOP_PER_THREAD" | bc))
        TIME_STEP_VALUE_NOW=$(echo "$PARALLEL_IT_ID * $NUM_PTS_IN_LOOP_PER_THREAD * $DT_VAL" | bc)
        TIME_RANGE_END_NOW=$(echo "($PARALLEL_IT_ID + 1) * $NUM_PTS_IN_LOOP_PER_THREAD * $DT_VAL" | bc)
        # Generate unique job name
        if [ "$NUM_PTS_IN_LOOP_PER_THREAD" -eq 1 ]; then
            JOBNAME="T${TIME_STEP_VALUE_NOW}_QD${MODE}.${SUBTYPE}_TMAX_${TMAXval}_Tunit_${T_UNIT}_Eq${quenched_E_val}_Jmax${JMAX}_Brot${BROT}_s${SPIN}_E${EL_FIELD}_B${MG_FIELD}_A${A_FACTOR}_del${D_EL}_dmg${D_MG}_Mtot${M_TOT}_NMOL${N_MOL}_${job_state_id}"
        else
            JOBNAME="T${TIME_STEP_VALUE_NOW}-${TIME_RANGE_END_NOW}_QD${MODE}.${SUBTYPE}_TMAX_${TMAXval}_Tunit_${T_UNIT}_Eq${quenched_E_val}_Rq_$(echo "$quenched_R_val" | tr ' ' '_')_Jmax${JMAX}_Brot${BROT}_s${SPIN}_E${EL_FIELD}_B${MG_FIELD}_A${A_FACTOR}_del${D_EL}_dmg${D_MG}_Mtot${M_TOT}_NMOL${N_MOL}_${job_state_id}"
        fi

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
            -e "s/BASIS_STATE_NO/${basis_state_no}/g" \
            -e "s/PRODUCT_STATE/${product_state}/g" \
            -e "s/IDSTR/${IDSTR}/g" \
            -e "s/DTNOW/${DTNOW}/g" \
            -e "s/ITSNUMOVERRIDE/${NUM_PTS_IN_LOOP_PER_THREAD}/g" \
            -e "s/ITNUMVAL/${IT_ID_NOW_INT}/g" \
            JuliaFewMoleculesQuenchDistributed > ${JOBNAME}.slurm

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
else
    if [ "$SUBTYPE" -eq 1 ]; then
        # Consolidation of distributed results
        JOBNAME="Consolidate_QD${MODE}.${SUBTYPE}_TMAX_${TMAXval}_Tunit_${T_UNIT}_Eq${quenched_E_val}_Jmax${JMAX}_Brot${BROT}_s${SPIN}_E${EL_FIELD}_B${MG_FIELD}_A${A_FACTOR}_del${D_EL}_dmg${D_MG}_Mtot${M_TOT}_NMOL${N_MOL}_${job_state_id}"

        echo "Creating SLURM job: $JOBNAME"

        sed -e "s/JOBNAME/${JOBNAME}/g" \
                -e "s/TIME/${TIME}/g" \
                -e "s/RAM/${NRAM_CONSOLIDATE}/g" \
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
                -e "s/ITSNUMOVERRIDE/${NUM_PTS_IN_LOOP_PER_THREAD}/g" \
                -e "s/ITNUMVAL/1/g" \
                JuliaFewMoleculesQuenchDistributed > ${JOBNAME}.slurm

        sbatch ${JOBNAME}.slurm

        if [ $? -eq 0 ]; then
                echo "✓ Job $JOBNAME submitted successfully"
            else
                echo "✗ Failed to submit job $JOBNAME"
            fi

        # Clean up job script
            rm -f ${JOBNAME}.slurm
    fi
fi



    

echo "========================================="
echo "All quench dynamics jobs submitted!"
echo "Finished at: " $(date)
echo "Total jobs submitted: $job_count"
echo "Monitor jobs with: squeue -u \$USER"
echo "Results will be in ${SCRATCH}/Julia_Few_Molecules/data/spectrum_results/ with organized hierarchy"
echo "========================================="
