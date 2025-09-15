# FewMolecules.jl

A high-performance Julia package for quantum mechanical calculations of few interacting magnetic polar molecules in external electric and magnetic fields. This package enables precise numerical studies of ultracold molecular physics, including energy spectra, quantum dynamics, and many-body interactions.

## Overview

FewMolecules.jl provides a comprehensive framework for modeling systems of 1-3 interacting magnetic polar molecules, in the description of Hund's case (e). It incorporates:

- **Full rotational and spin degrees of freedom** with arbitrary J (rotational) and S (total electronic) quantum numbers
- **Realistic molecular interactions**: electric and magnetic dipole-dipole interactions, spin-rotation coupling
- **External field effects**: arbitrary electric and magnetic field orientations and strengths
- **Time evolution dynamics**: quantum quench simulations and real-time evolution
- **High-performance computing**: multi-threaded calculations with optimized linear algebra
- **Flexible geometry**: arbitrary molecular arrangements in 2D/3D space

The package is designed for studies of ultracold molecular gases, quantum magnetism and few-body physics.

## Key Features

### Physical Models
- **Molecular Hamiltonians**: Complete treatment of magnetic polar molecules with spin-rotation coupling
- **Intermolecular Interactions**: 
  - Dipole-dipole interactions: Electric and Magnetic (including anisotropic terms)
  - Customizable interaction strengths and selection rules
- **External Fields**: 
  - Arbitrary electric field strength and orientation
  - Arbitrary magnetic field strength and orientation  
  - Field-dependent molecular polarization

    **Current limitation:** Electric and Magnetic fields must be parallel.

- **Quantum Numbers**: Support for arbitrary rotational angular momentum $j$ and total electronic spin $s$ quantum numbers.

    **Current limitation:** For now all molecules are assumed to have the same maximal quantum numbers.

### Computational Capabilities
- **Spectrum Calculations**: Energy eigenvalues and eigenstates as functions of external parameters
- **Time Evolution**: Real-time quantum dynamics with various time units (harmonic, SI units)
- **Quantum Quenches**: Sudden parameter changes and non-equilibrium dynamics
- **Observable Calculations**: Expectation values of individual molecular angular momentum operators and whole system observables
- **Parameter Sweeps**: Automated calculations over ranges of physical parameters

### Performance Features
- **Multi-threading**: Efficient parallelization using Julia's threading capabilities
- **Distributed calculation capabilities:** The code was fitted to allow for distributed parallel execution on a computer cluster as batches of independent calculation jobs. 
- **Memory Optimization**: Sparse matrix representations and efficient storage
- **Numerical Precision**: Configurable floating-point precision (Float32, Float64, ComplexF32, ComplexF64)
- **Progress Monitoring**: Real-time progress bars and status reporting
- **Profiling Support**: Built-in performance profiling and benchmarking

## Scientific Applications

This package has been used to study:
- **Few-body physics**: one-, two- and three-molecule energy spectra and formation of anticrossings under changes of simulation parameters
- **Non-equilibrium dynamics**: Quantum quenches of states undergoing anticrossings.

## Installation

### Prerequisites
- **Julia**: Version 1.8+ recommended
- **BLAS/LAPACK**: Optimized linear algebra libraries (automatically handled by Julia)

### Setup
1. Clone the repository:
   ```bash
   git clone https://github.com/your-org/Few-Magnetic-Polar-Molecules.git
   cd Few-Magnetic-Polar-Molecules
   ```

2. Install Julia dependencies:
   ```bash
   julia --project=. -e "using Pkg; Pkg.instantiate()"
   ```

3. (Optional) Verify installation:
   ```bash
   ./local_execution_scripts/julia_test_local.sh
   ```

## Quick Start

### Basic Spectrum Calculation
Calculate energy spectrum vs electric field in the range $\mathcal{E} \in [0 - 15]\,\mathrm{kV/cm}$ for two molecules:

```bash
julia --project=. src/FewMolecules.jl 2.0 "100 0" --B 50 --E 0 --d_el 0.1 --d_mg 10 --A 100 --Mtot -4.0 --START_IT 1 --END_IT 101 --step 0.15 --looped-eigenproblem
```

### Time Evolution Simulation
Run quantum quench dynamics:

```bash
julia --project=. src/FewMolecules.jl 1.0 "100 0" --B 50 --quenched_E 11.49914 --productstate "1,-1,2.5,-2.5,1,-1,2.5,0.5" --tmax 1000 --dt 0.25
```

### Using Convenience Scripts
Pre-configured scripts for common calculations:

```bash
# Quick functionality test
./scripts/julia_test_local.sh

# Spectrum calculation with parameter sweeps  
./scripts/julia_spectrum_local.sh

# Quantum quench dynamics
./scripts/julia_quench_local.sh
```

For further examples of spectrum and quench dynamics calculation scripts, please have a look into the [local_execution_scripts](./local_execution_scripts/) and [slurm_scripts](./slurm_scripts/) folders. Example SLURM script configured for PLGrid Ares supercomputer that showcase the distributed calculation capabilities can be found in [slurm_scripts](./slurm_scripts/) directory.

## Usage Guide

### Command Line Interface

The main program accepts the following syntax:
```bash
julia --project=. src/FewMolecules.jl MODE.SUBTYPE GEOMETRY [OPTIONS...]
```

#### Modes and Subtypes
- **Mode 0**: Eigenenergies calculation

    Allows for use of optimized routines for hermitian matrices, reducing computational cost. Recommended for large parameter sweeps, where eigenstates are not essential.

- **Mode 2**: Parameter sweep calculations

    Full eigenproblem (eigenvalues + eigenstates) calculation. The LAPACK hermitian optimized routines for hermitian matrices treatment are currently not implemented in Julia, so it is more computationally intensive than pure eigenvalue calculation. Recommended for calculation of the regions of interest that will benefit from the eigenstate breakdown for detailed analysis.

    **The two modes share common subtypes:**

  - Subtype 0: Spectrum vs electric field $\mathcal{E}$,
  - Subtype 1: Spectrum vs magnetic field $\mathcal{B}$,
  - Subtype 2: Spectrum vs A-factor (spin-rotation coupling $\gamma$ ),
  - Subtype 3: Spectrum vs distance $R_{ij}$ scaling,
  - Subtype 4: Spectrum vs electric dipole moment $d_o$,
  - Subtype 5: Spectrum vs magnetic dipole moment $\mu_B$,
- **Mode 1**: Time evolution dynamics  

    Time dynamics under rapidly changed (quenched) conditions of electric, magnetic fields and / or system geometry. 

  - Subtype 0: Calculation of the time dynamics of an initial state under quenched conditions.
    
    The following two subtypes are implemented but not incorporated into the package yet:

  - Subtype 1: consolidate results from distributed calculation of quench dynamics on a computer cluster,
  - Subtype 2: plot result time series and Fourier transforms.

#### Geometry Specification
Molecular positions are specified as space-separated values:
- **1 molecule**: (implicit at origin)
- **2 molecules**: `"R01 Φ01"` (R in nm, Φ in degrees)  
- **3 molecules**: `"R01 Φ01 R02 Φ02"` 

Examples:
- `"200 90"` - Two molecules: one at origin, one at 200 nm at 90°
- `"500 30 500 90"` - Three molecules in triangular arrangement

### Key Parameters

#### Physical Parameters
- `--E FIELD` - Electric field strength (kV/cm)
- `--B FIELD` - Magnetic field strength (Gauss)
- `--d_el DIPOLE` - Electric dipole moment (Debye)
- `--d_mg DIPOLE` - Magnetic dipole moment (μB)
- `--A FACTOR` - Spin-rotation coupling constant (kHz)
- `--Mtot PROJECTION` - Total angular momentum projection (float)

#### Field Orientations
- `--theta_el ANGLE` - Electric field polar angle $\Theta_\mathcal{E}$ (degrees)  
- `--phi_el ANGLE` - Electric field azimuthal angle $\Phi_\mathcal{E}$ (degrees)
- `--theta_mg ANGLE` - Magnetic field polar angle $\Theta_\mathcal{B}$ (degrees)
- `--phi_mg ANGLE` - Magnetic field azimuthal angle $\Phi_\mathcal{B}$ (degrees)

#### Quantum Numbers
- `--JMAX J` - Maximum total angular momentum for each particle $j$ (int)
- `--spin S` - Total electronic angular momentum quantum number ($s \equiv s/2$, so --spin 13 $\equiv s=13/2$)
- `--Brot CONSTANT` - Rotational constant $B$ (GHz)

#### Numerical Settings
- `--hamiltonian-dtype TYPE` - Numerical precision (Float64, ComplexF64)
- `--basis-method METHOD` - Basis construction method
- `--hamiltonian-method METHOD` - Hamiltonian construction method
- `--eigensolver SOLVER` - Eigenvalue solver selection

#### Time Evolution Parameters  
- `--tmax TIME` - Maximum evolution time
- `--dt TIMESTEP` - Time step size
- `--time-unit UNIT` - Time units (harmonic, s, ms, μs, ns, ps, fs)
- `--quenched_E FIELD` - Post-quench electric field
- `--quenched_R GEOMETRY` - Post-quench geometry
- `--productstate STATE` - Initial product state specification

## Architecture

### Core Components

#### `FewMolecules.jl`
- The main entry point of the package. Handles command-line argument parsing, initialization, and dispatches to the appropriate calculation mode (spectrum, quench, parameter sweeps, etc.).
- Loads all other modules and manages global configuration, BLAS backend selection, and profiling.
- Exports the main function for direct execution via `julia src/FewMolecules.jl ...` or module loading with `using FewMolecules` construct in Julia file or REPL.

#### `basis_cpp_faithful.jl`
- Implements the basis state generation for the system, faithfully ported from the original C++ algorithm.
- Handles quantum number conventions, total angular momentum projection constraints, and supports 1–4 molecules.
- Provides conversion utilities between C++-style (integer, s=2s) and Julia-style (float s=s) basis representations .

#### `eigenproblem.jl`
- Contains routines for solving the eigenvalue problem for the system Hamiltonian.
- Supports in-place diagonalization, unit conversions for eigenvalues, and phase standardization of eigenvectors.
- Used for both spectrum calculations and as a backend for time evolution.

#### `looped_eigenproblem.jl`
- Provides high-level routines for sweeping a parameter (electric/magnetic field, spin-rotation, geometry, dipole moments) and efficiently reusing the base Hamiltonian.
- Handles the logic for iterating over parameter ranges, updating only the relevant Hamiltonian terms, and saving results.
- Integrates progress bars, memory profiling, and result organization.

#### `hamiltonian_cpp_faithful.jl`
- Constructs the full system Hamiltonian using a modular approach, closely following the original C++ implementation.
- Adds rotational, electric field, magnetic field, spin-rotation, electric dipolar, and magnetic dipolar terms as needed.
- Handles geometry and field orientation, supports sparse matrix storage, and ensures Hermiticity.

#### `ham_functions_cpp_faithful.jl`
- Contains helper functions for Hamiltonian construction, including Clebsch-Gordan coefficients, dipole and spin operators, and selection rules.
- Implements the detailed logic for matrix elements, selection rules, and Hermitization for both dense and sparse matrices.

#### `parameters.jl`
- Defines the `Params` struct, which holds all simulation parameters, including geometry, fields, molecular properties, and numerical settings.
- Handles parsing and validation of command-line arguments, unitful quantities, and parameter consistency.
- Provides helpers for parameter iteration and skip-term logic.

#### `PhysicalConstantsJL.jl`
- Provides a clear, centralized definition of all physical constants used in the package, using the CODATA 2022 values.
- Exports constants such as the fine-structure constant, Bohr radius, Planck constant, and more, all as Unitful quantities.

#### `units.jl`
- Implements all unit conversion logic, including conversions to and from atomic units, or SI.
- Provides functions for converting frequencies, wavenumbers, dipole moments, fields, and times to atomic units.
- Ensures consistent handling of units throughout the codebase.

#### `quench.jl`
- Implements the time evolution (quench dynamics) logic, including initial state preparation, time propagation, and observable calculation.
- Supports both product state, eigenstate and basis state initializations, and handles arbitrary observables.
- Includes parallelized and vectorized routines for efficient time evolution and result saving.

#### `profiling.jl`
- Provides lightweight wrappers for memory and performance profiling.
- Tracks memory allocations, peak usage, and timing for key operations, with reporting utilities for performance analysis.
- Can be enabled or disabled as needed for debugging or optimization.

#### `utility.jl`
- Contains general-purpose helper functions for file naming, directory management, parameter formatting, and type parsing.
- Handles the construction of hierarchical data directories and filenames based on simulation parameters.
- Provides helpers for formatting, parameter string conversion, and result organization.

#### `saving_utils.jl`
- Implements routines for saving eigenvalues, eigenvectors, simulation parameters, and basis states in various formats (Apache Arrow, JLD2, text).
- Handles thresholding, sparse storage, and metadata for efficient and organized result storage.
- Ensures compatibility with downstream analysis tools.

## Output and Results

### Data Organization
Results are automatically organized in the `data/` directory:
```
data/
├── spectrum_results/    # Energy spectra and parameter sweeps
├── quench_results/      # Time evolution data
└── logs/               # Calculation logs and metadata
```

### File Formats

- **Eigenvalues and spectra**: Saved in [Apache Arrow](https://arrow.apache.org/) format (`.arrow`) for efficient columnar storage and fast loading, or optionally as plain text (`.txt`).
- **Eigenstates (eigenvectors)**: Saved as Arrow files (dense encoding), or as [JLD2](https://github.com/JuliaIO/JLD2.jl) (`.jld2`) files for large or sparse data. The `JLD2` format SparseVectors saving allows to save 10-100x the memory compared to saving dense eigenvectors in `.arrow` or `.csv` format. Text output is also supported for human-readable inspection.
- **Simulation parameters and basis**: Saved as [JLD2](https://github.com/JuliaIO/JLD2.jl) (`.jld2`) or Arrow files.
- **Quench dynamics and time evolution**: Saved as plain text or `.csv` files, with headers describing the observable and simulation parameters.
- **Metadata**: Parameters and run information are saved alongside results in JLD2 or Arrow files for reproducibility.

### Visualization
While the package focuses on computation, output formats are compatible with:
- **Julia Plotting**: Plots.jl, PlotlyJS.jl for interactive visualization
- **Python**: Easy import into matplotlib, numpy for analysis
- **Mathematica**: Direct import of numerical results
- **MATLAB**: Compatible data formats for further analysis

## Performance Considerations

### Scaling
- **System Size**: Computational cost scales exponentially with number of molecules
- **Basis Size**: Scales as (2J+1)^N × (2S+1)^N for N molecules
- **Memory Usage**: Dominated by Hamiltonian matrix storage O(basis_size^2)
- **Time Evolution**: Cost per time step scales as O(basis_size^2)

### Optimization Tips
- Use appropriate `JULIA_NUM_THREADS` for your system
- Choose minimal `JMAX` and `SPIN` values for your problem
- Consider `Float32` precision for large calculations if accuracy permits
- Use sparse representations for systems with good quantum number conservation

### Hardware Requirements
- **RAM**: 2GB minimum for 1-molecule calculations, 10GB minimum for 2-molecule calculations and 32GB+ recommended for 3-molecule calculations
- **CPU**: Multi-core processor. A substantial performance increase is visible if the processor is compatible with [Intel® oneAPI Math Kernel Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) or [Apple Accelerate](https://developer.apple.com/accelerate/) architecture.
- **Storage**: Fast SSD recommended for I/O intensive parameter sweeps


<!-- ## Contributing

We welcome contributions to FewMolecules.jl! Please see our contributing guidelines:

1. **Issues**: Report bugs and request features via GitHub issues
2. **Pull Requests**: Submit improvements with tests and documentation
3. **Documentation**: Help improve examples and tutorials
4. **Benchmarking**: Contribute performance benchmarks and optimizations

### Development Setup
```bash
git clone https://github.com/your-org/Few-Magnetic-Polar-Molecules.git
cd Few-Magnetic-Polar-Molecules
julia --project=. -e "using Pkg; Pkg.develop()"
``` -->

## Citation

If you use FewMolecules.jl in your research, please cite:

```bibtex
@software{fewmolecules_jl,
  title={FewMolecules.jl: A Julia Package for Few-Body Quantum Calculations of Magnetic Polar Molecules},
  author={[Kacper Cybiński, Anna Dawid, Michał Suchorowski, Michał Tomza]},
  year={2024},
  url={https://github.com/your-org/Few-Magnetic-Polar-Molecules},
  doi={[DOI to be added after Zenodo submission]}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- **Documentation**: Check this README and inline documentation
- **Issues**: Submit bug reports and feature requests via GitHub issues  
- **Discussions**: Use GitHub discussions for questions and community support
- **Contact**: [kacper.cybinski@fuw.edu.pl](malito:kacper.cybinski@fuw.edu.pl) for direct inquiries

## Acknowledgments

This project makes extensive use of the Julia scientific computing ecosystem. In particular, we acknowledge:

- The [Julia language](https://julialang.org/) and its core developers ([MIT License](https://github.com/JuliaLang/julia/blob/master/LICENSE.md)).  
  If you use Julia in your research, please cite it as described in [Julia's citation guidelines](https://julialang.org/citation/).
- [LinearAlgebra.jl](https://github.com/JuliaLang/LinearAlgebra.jl) and [SparseArrays.jl](https://github.com/JuliaSparse/SparseArrays.jl) from the Julia standard library for high-performance dense and sparse matrix operations.
- [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) and [UnitfulAtomic.jl](https://github.com/sostock/UnitfulAtomic.jl) for physical units and atomic units support.
- [PhysicalConstants.jl](https://github.com/JuliaPhysics/PhysicalConstants.jl) for CODATA 2022 physical constants.
- [ArgParse.jl](https://github.com/carlobaldassi/ArgParse.jl) for command-line argument parsing.
- [Arrow.jl](https://github.com/JuliaData/Arrow.jl) and [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) for efficient data serialization.
- [ProgressBars.jl](https://github.com/JuliaCollections/ProgressBars.jl) for progress reporting.
- [WignerSymbols.jl](https://github.com/Jutho/WignerSymbols.jl) for angular momentum algebra.
- [MKL.jl](https://github.com/JuliaLinearAlgebra/MKL.jl) and [AppleAccelerate.jl](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl) for optimized BLAS backends (where available).

We thank the developers and maintainers of these packages for their contributions to the Julia ecosystem.  
If you use this software in your research, please also cite Julia and the relevant packages as described in their documentation.

**Special thanks to collaborators who provided feedback and testing during development.** <br>
This goes especially to Anna Dawid (LIACS @ Leiden University), who created the original C++ implementation of the core functionalities of this package. Further huge thanks are to Michał Suchorowski (University of Warsaw) for continued support and endless hours of discussions throughout the development process, as well as Michał Tomza (University of Warsaw) for providing his expertise and for the oversight of the project.

---

*FewMolecules.jl - Enabling precision quantum calculations for ultracold molecular physics*
