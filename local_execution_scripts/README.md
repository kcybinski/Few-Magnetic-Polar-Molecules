# FewMolecules.jl Scripts

This directory contains shell scripts for simplified execution of the FewMolecules.jl package. These scripts provide convenient wrappers around the main Julia program with pre-configured parameters for common calculation types.

## Overview

The scripts are designed to streamline the execution of FewMolecules.jl calculations by providing:
- Pre-configured parameter sets for different calculation modes
- Automatic directory setup and navigation
- Julia multithreading configuration
- Progress reporting and error handling

**Note**: These scripts are not the exclusive way to run FewMolecules.jl. You can always run the Julia package directly with custom parameters.

## Available Scripts

### 1. `julia_test_local.sh`
**Purpose**: Quick verification that the FewMolecules.jl package is working correctly.

**What it does**:
- Runs a simple spectrum calculation (Mode 0, Subtype 0)
- Uses minimal computational parameters for fast execution
- Validates Julia installation and project structure
- Provides basic functionality testing

**Use when**: Setting up the environment or troubleshooting installation issues.

### 2. `julia_spectrum_local.sh`
**Purpose**: Perform spectrum calculations across parameter ranges.

**What it does**:
- Executes spectrum calculations (Mode 2) with various subtypes
- Supports sweeping over electric field, magnetic field, A-factor, geometry scaling, and dipole moments
- Configured for production-level calculations with fine parameter steps
- Includes comprehensive parameter documentation

**Use when**: Computing energy spectra as functions of external parameters.

### 3. `julia_quench_local.sh`
**Purpose**: Simulate quench dynamics of molecular systems.

**What it does**:
- Runs time-evolution calculations (Mode 1) after sudden parameter changes
- Supports both harmonic and SI time units
- Calculates observables (M, Ms, m, ms) during time evolution
- Includes extensive parameter sets for initial and quenched states

**Use when**: Studying non-equilibrium dynamics and quantum quenches.

## Parameter Customization

### Important Notes:
- **All parameters in the scripts are for demonstration purposes only**
- **Modify parameters according to your specific research needs**
- **The current values are examples and may not represent physical systems of interest**

### Key Parameters to Customize:

#### Physical Parameters:
```bash
D_EL='0.1'        # Electric dipole moment (Debye)
D_MG='10'         # Magnetic dipole moment (μB)  
A_FACTOR='100'    # Spin-rotation coupling (kHz)
MG_FIELD='50'     # Magnetic field (Gauss)
EL_FIELD='11.495' # Electric field (kV/cm)
```

#### Geometry Configuration:
```bash
GEOMETRY='200 90'              # Two molecules: R₀₁=200nm, Φ₀₁=90°
# GEOMETRY='500 30 500 90'     # Three molecules (uncomment to use)
```

#### Quantum Numbers:
```bash
JMAX='2'          # Maximum total angular momentum
SPIN='5'          # Spin quantum number  
MTOT='-4.0'       # Total angular momentum projection
```

#### Computational Settings:
```bash
JULIA_NUM_THREADS=11          # CPU threads for distributed execution
DTYPE="Float64"               # Numerical precision
```

## Usage Instructions

1. **Make scripts executable**:
   ```bash
   chmod +x *.sh
   ```

2. **Run a script**:
   ```bash
   ./julia_test_local.sh        # Quick test
   ./julia_spectrum_local.sh    # Spectrum calculation
   ./julia_quench_local.sh      # Quench dynamics
   ```

3. **Monitor progress**:
   - All scripts include progress bars and status reporting
   - Results are saved to `data/` subdirectories
   - Check console output for completion status

## Customization Guide

### To modify parameters:
1. Open the desired script in a text editor
2. Locate the parameter section (clearly marked with comments)
3. Modify values according to your research needs
4. Save and run the script

### To add new parameter sweeps:
- Modify `START_IT`, `END_IT`, and `STEP` variables
- Adjust the subtype to sweep different parameters
- See FewMolecules.jl documentation for available options

### To change calculation modes:
- Set `MODE` and `SUBTYPE` variables appropriately
- Refer to the main package documentation for mode descriptions

## Output

Results are typically saved to:
- `data/spectrum_results/` (spectrum calculations)
- `data/quench_results/` (quench dynamics)
- Individual files are named based on calculation parameters

## Troubleshooting

- **Julia not found**: Ensure Julia is installed and in your PATH
- **Package not found**: Run `julia --project=. -e "using Pkg; Pkg.instantiate()"`
- **Memory issues**: Reduce `JULIA_NUM_THREADS` or calculation parameters
- **Slow performance**: Check that multithreading is working (`JULIA_NUM_THREADS`)

## Direct Julia Usage

These scripts are convenience wrappers. You can always run FewMolecules.jl directly:

```bash
julia --project=. src/FewMolecules.jl [MODE].[SUBTYPE] [GEOMETRY] [OPTIONS...]
```

See the main package documentation for complete parameter lists and advanced options.
