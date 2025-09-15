module FewMolecules

using ArgParse
using Dates
using LinearAlgebra # Added for BLAS.vendor() in __init__
using Unitful # For unit handling, if needed

include("parameters.jl")                    # Defines module Parameters
include("PhysicalConstantsJL.jl")           # Defines module PhysicalConstantsJL
include("units.jl")                         # Defines module Units
include("utility.jl")                       # Defines module Utility
include("profiling.jl")                     # Defines module Profiling
include("saving_utils.jl")                  # Defines module SavingUtils
include("basis_cpp_faithful.jl")            # Defines module BasisCppFaithful
include("ham_functions_cpp_faithful.jl")    # Needed for HamiltonianCpp
include("hamiltonian_cpp_faithful.jl")      # Defines module HamiltonianCpp
include("eigenproblem.jl")                  # Defines module Eigenproblem
include("looped_eigenproblem.jl")           # Defines module LoopedEigenproblem
include("quench.jl")                        # Defines module Quench

# Export main to be callable
export main
export MKL_AVAILABLE
export current_params # Export current_params for access in main

# Global variable to track MKL availability
const MKL_AVAILABLE = Ref(false)
function __init__()
    if Sys.isapple()
        try
            @eval using AppleAccelerate
            println("INFO: Running on macOS. Attempting to use AppleAccelerate for BLAS operations.")
            # AppleAccelerate typically sets itself as the BLAS provider automatically upon loading.
            println("INFO: BLAS config after attempting to load AppleAccelerate: ", LinearAlgebra.BLAS.get_config())
            MKL_AVAILABLE[] = false  # AppleAccelerate is not MKL
        catch e
            println("WARNING: AppleAccelerate.jl is in dependencies but could not be loaded. Error: ", e)
            println("INFO: Continuing with default BLAS on macOS: ", LinearAlgebra.BLAS.get_config())
            MKL_AVAILABLE[] = false
        end
    else # For non-macOS systems (like Linux on Ares)
        try
            @eval using MKL
            # MKL.set_mkl_threads(1) # Set MKL to use single thread for consistency
            println("INFO: Running on Linux. Attempting to use Intel MKL for BLAS operations.")
            println("INFO: BLAS config after MKL.set_mkl_threads(1): ", LinearAlgebra.BLAS.get_config())
            MKL_AVAILABLE[] = true  # MKL successfully loaded
        catch e
            println("WARNING: MKL.jl is in dependencies but could not be loaded or MKL.set_mkl_threads(1) failed. Error: ", e)
            println("INFO: Continuing with default BLAS on Linux system: ", LinearAlgebra.BLAS.get_config())
            MKL_AVAILABLE[] = false
        end
    end
end

function parse_cli_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "mode"
            help = """
            Calculation mode: 
            0.X (spectrum), 1.X (quench), 2.X (spectrum + eigenstates), 3.X (magnetization), 4.X (print eigenfunctions)
            Subtypes:
            [0.X / 2.X] Spectrum subtype: 0-vs Electric field, 1-vs Magnetic field, 2-vs A factor, 3-vs Electric dipole, 4-vs Magnetic dipole
            [1.X] Quench subtype: 0-vs calculate quench , implementation started, but not complete: 1-vs consolidate distributed results, 2-vs plot result time series and Fourier transforms
            """
            arg_type = Float64
            required = true
        "geometry..."
            help = "Geometry parameters: -R_01 -Φ_01 -R_02 -Φ_02 -R_03 -Φ_03 (max 4 molecules, leave empty for single molecule)"
            nargs = '*' # Allow zero or more geometry parameters
            arg_type = Float64
        "--E"
            arg_type = Float64
            default = 0.0
            help = "Electric field value (e.g., 10.0 for 10 kV/cm)"
        "--E-unit"
            arg_type = String
            default = "kV/cm"
            help = "Electric field unit (e.g., kV/cm, V/cm)"
        "--phi_el"
            arg_type = Float64
            default = 0.0
            help = "Phi_E (deg)"
        "--theta_el"
            arg_type = Float64
            default = 0.0
            help = "Theta_E (deg)"
        "--B"
            arg_type = Float64
            default = 0.0
            help = "Magnetic field (G)"
        # No --B-unit, assuming always Gauss as per units.jl
        "--phi_mg"
            arg_type = Float64
            default = 0.0
            help = "Phi_B (deg)"
        "--theta_mg"
            arg_type = Float64
            default = 0.0
            help = "Theta_B (deg)"
        "--A"
            arg_type = Float64
            default = 0.0
            help = "Spin-rotation coupling value (e.g., 50.0 for 50 kHz)"
        "--A-unit"
            arg_type = String
            default = "kHz"
            help = "Spin-rotation coupling unit (e.g., kHz, MHz, Hz)"
        "--d_el"
            arg_type = Float64
            default = 0.0
            help = "d_el (D)"
        "--d_mg"
            arg_type = Float64
            default = 0.0
            help = "d_mg (mu_B)"
        "--Mtot"
            arg_type = Float64 # Changed from Int to Float64
            default = 0.0
            help = "Mtot (physical total projection, can be half-integer)"
        "--START_IT"
            arg_type = Int
            default = 0
            help = "START_IT"
        "--END_IT"
            arg_type = Int
            default = 0
            help = "END_IT"
        "--step"
            arg_type = Float64
            default = 1
            help = "Looped Eigenproblem step size"
        "--IT_NO_OVERRIDE"
            arg_type = Int
            default = nothing
            help = "Perform a single step from the looped calculation (added for cluster parallel computation capability)."
        "--ITS_PER_OVERRIDE"
            arg_type = Int
            default = nothing
            help = "Extends the multithreaded loop to N iterations (added for cluster parallel computation capability)."
        "--JMAX"
            arg_type = Int
            default = 2 
            help = "JMAX"
        "--QN_per_mol"
            arg_type = Int
            default = 4  # Default: 4 quantum numbers per molecule (j, m, S, ms)
            help = "Number of quantum numbers per molecule"
        "--Brot"
            arg_type = Float64
            default = 0.29  # C++ default: B_rot = 0.29 GHz (KRb: 0.5264, DyEr: 0.29)
            help = "B_rot value (e.g., 20.0 for 20 GHz)"
        "--Brot-unit"
            arg_type = String
            default = "GHz"
            help = "B_rot unit (e.g., GHz, MHz, Hz)"
        "--spin"
            arg_type = Float64
            default = 13.0  # C++ default: S1 = 13 (*2: ErA = 13/2, DyEr = 17/2)
            help = "Spin quantum number (input 2*S, e.g., 1 for S=1/2, 2 for S=1)"
        "--productstate"
            help = "Initial product state (comma-separated)"
        "--quenched_E"
            arg_type = Float64
            help = "Quenched E (kV/cm)"
        "--quenched_B"
            arg_type = Float64
            help = "Quenched B (G)"
        "--quenched_R"
            arg_type = Float64
            nargs = '+'
            default = nothing
            help = "Quenched R (nm)"
        "--quench-observables"
            nargs = '+'
            arg_type = String
            help = """
List of observables to track during quench dynamics. Possible options:
  Ms  : total spin projection (sum of ms1, ms2, ...)
  M   : total orbital projection (sum of m1, m2, ...)
  ms1 : spin projection for molecule 1
  ms2 : spin projection for molecule 2
  ms3 : spin projection for molecule 3
  m1  : orbital projection for molecule 1
  m2  : orbital projection for molecule 2
  m3  : orbital projection for molecule 3
Any combination is allowed, e.g. --quench-observables Ms M ms1 m2
Assuming QNs per molecule = 4, for 2 molecules: [j1 m1 S1 ms1 j2 m2 S2 ms2]
"""
        "--harmonic-frequency"
            arg_type = Float64
            default = 100.0
            help = "Trap harmonic frequency value for time units (e.g., 1000.0 for 1000 MHz), to be multiplied by 2π. The default is 100 kHz"
        "--harmonic-frequency-unit"
            arg_type = String
            default = "kHz"
            help = "Unit for harmonic frequency (Hz, kHz, MHz, GHz, etc.). The default is kHz."
        "--harmonic-frequency-source"
            arg_type = String
            default = "Brot"
            help = "Source of harmonic frequency: 'Brot' (use rotational constant), 'A' (use spin-rotation), 'custom' (use --harmonic-frequency)"
        "--time-unit"
            arg_type = String
            default = "harmonic"
            help = "Input time unit: 'harmonic' (dimensionless 2π/ω units) or SI units ('s', 'ms', 'μs', 'ns', 'ps', 'fs')"
        "--tmax"
            arg_type = Float64
            default = 1.0
            help = "Maximum evolution time (unit specified by --time-unit)"
        "--dt"
            arg_type = Float64
            default = 0.01
            help = "Time step (unit specified by --time-unit)"
        "--outputdir"
            help = "Output directory"
        "--basis-method"
            arg_type = String
            default = "cpp-faithful"
            help = "Basis generation method: 'julia' (deprecated, currently not supported) or 'cpp-faithful' (default, C++ faithful port for validation)"
        "--hamiltonian-method"
            arg_type = String
            default = "cpp-faithful"
            help = "Hamiltonian construction method: 'julia' (Julia implementation - deprecated, currently not supported) or 'cpp-faithful' (C++ faithful port for validation)"
        "--save-hamiltonian"
            action = :store_true
            help = "Save Hamiltonian matrix to file for comparison with C++ implementation"
        "--hamiltonian-dtype"
            arg_type = String
            default = "ComplexF64" # Default to ComplexF64 for Hamiltonian matrix
            help = "Data type for Hamiltonian matrix: 'ComplexF64' (default), 'Float64' (real-valued Hamiltonian) or reduced 'Float32' or 'ComplexF16' (for smaller memory footprint)"
        "--skip-terms"
            arg_type = String
            default = ""
            help = "Comma-separated list of terms to artificially skip in Hamiltonian construction. Available options are: rotation, electric_field, magnetic_field, spin_rotation, electric_dipolar, magnetic_dipolar."
        "--eigensolver"
            arg_type = String
            default = "auto"
            help = "Eigensolver method: 'auto' (automatic selection), 'mkl' (force MKL), 'julia' (force Julia built-in)"
        "--spectra-unit"
            arg_type = String
            default = "Brot"
            help = "Unit for eigenvalue output: 'Brot' (B_rot units, default), 'Hz' (frequency), 'Hartree' (atomic units)"
        "--looped-eigenproblem"
            action = :store_true
            help = "Switches the task 0 (diagonalization) from single realization to looped implementation between START_IT and END_IT with step."
        "--profile"
            action = :store_true
            help = "Enable basic performance and memory profiling"
        "--initial-eigenstate"
            arg_type = Int
            default = 1
            help = "Index of initial eigenstate to use for quench (1 = ground state, 2 = first excited, ...)"
        "--initial-basis-state"
            arg_type = Int
            default = nothing
            help = "Index of initial basis state to use for quench (1 = first basis state, 2 = second basis state, ...)"
        "--colored-print"
            action = :store_true
            help = "Enable colored prints in terminal (not advidsed for cluster calculations, print to .log might get corrupted)."
        "--progressbar"
            action = :store_true
            help = "Toggles interim printing in favour of a progress bar in looped eigenproblem calculation."
        "--run-id-string"
            arg_type = String
            help = "When running this from a shell script, pass to this argument the date/time of run start."
        "--save-root-folder"
            arg_type = String
            default = "./"
            help = "Root folder for saving results. Default is current directory './'. On a cluster, pass the scratch directory."
    end
    return parse_args(s) # ArgParse uses ARGS by default
end

function main()
    args = parse_cli_args()

    # Initialize simple profiling
    profile_enabled = get(args, "profile", false)
    Profiling.init_profiling(profile_enabled)
    
    Profiling.@time_operation "Parameter_Initialization" Parameters.init_from_args(args)
    # Access mode from parameters after initialization
    current_params = Parameters.PARAMS[]
    mode = current_params.mode
    vsWhat = current_params.vsWhat
    
    current_params.profile = profile_enabled
    
    # Track parameters memory
    Profiling.track_memory("Parameters", Profiling.get_memory_size(Parameters.PARAMS[]))
    
    flush(stdout)
    if mode == 0 || mode == 2
        SavingUtils.save_simulation_params(current_params)

        println("Running spectrum calculation...")
        # Save simulation parameters to file
        
        # Basis generation with performance tracking
        basis_method = args["basis-method"]
        hamiltonian_method = args["hamiltonian-method"]
        
        channels = Profiling.@time_operation "Basis_Generation" begin
            if basis_method == "cpp-faithful"
                println("INFO: Using C++ faithful basis generation (validation/debugging)")
                BasisCppFaithful.generate_channels_cpp_style(
                    current_params.Mtot,
                    false, # print_debug
                    JMAX=current_params.JMAX,
                    spin=current_params.spin / 2.0,  # Convert 2*S to S
                    geometry=args["geometry..."] # Pass geometry from args
                )
            else
                error("Invalid basis method: '$basis_method'. Valid option is 'cpp-faithful'")
            end
        end

        mem_req(s, ::Type{T}) where {T} = 4s^2*sizeof(T) / (1 << 30)

        T = Utility.parse_hamiltonian_dtype(args["hamiltonian-dtype"])

        println("INFO: Estimated memory required for dense Hamiltonian matrix of size $(length(channels))×$(length(channels)) with type $(T): $(mem_req(length(channels), T)) GB")

        SavingUtils.save_basis(channels, current_params)

        QN = length(channels)
        println("INFO: Number of channels generated: $QN")
        
        flush(stdout)
        # Track channels memory
        Profiling.track_memory("Channels", Profiling.get_memory_size(channels))
        
        if QN == 0
            println("Warning: No channels generated. Hamiltonian will be empty. Check JMAX, spin, Mtot, and geometry.")
            # Generate clean parameter strings for consistent filename format
            filename = Utility.create_filename(current_params, 0, "eigenvalues_empty")
            open(filename, "w") do io
                println(io, "# No channels generated for these parameters.")
            end
            println("Spectrum calculation skipped. Empty results file created: $filename")
        else
            # Check if looped eigenproblem is requested
            if args["looped-eigenproblem"]
                println("Running looped eigenproblem calculation...")
                
                # Determine how many eigenstates to save
                how_many_eigenstates_to_save = QN  # Default to basis size
                
                # Create looped parameters with additional solver settings
                looped_params = deepcopy(current_params)
                
                flush(stdout)
                # Run looped eigenproblem with direct parameter passing
                Profiling.@time_operation "Looped_Eigenproblem" begin
                    LoopedEigenproblem.looped_eigenproblem(
                        channels, QN, looped_params, how_many_eigenstates_to_save;
                        param_strings=Utility.params_to_cpp_strings(current_params),
                        spectra_unit=args["spectra-unit"],
                        eigensolver_method=args["eigensolver"],
                        mkl_available=MKL_AVAILABLE[]
                    )
                end
            else
                # Hamiltonian construction with detailed step tracking
                H_sparse = Profiling.@time_operation "Hamiltonian_Construction_Total" begin
                    if hamiltonian_method == "cpp-faithful"
                        println("INFO: Using C++ faithful Hamiltonian construction (production code)")
                        HamiltonianCpp.buildHamiltonian(channels, QN, current_params)
                    else
                        error("Invalid hamiltonian method: '$hamiltonian_method'. Valid option is 'cpp-faithful'")
                    end
                end
                # Run single diagonalization
                println("Running single diagonalization...")
                
                # Convert sparse to dense for single diagonalization
                H = Matrix(H_sparse)

                # Hermiticity check with performance tracking
                Profiling.@time_operation "Hermiticity_Check" begin
                    if !ishermitian(H)
                        println("Warning: Hamiltonian matrix is not Hermitian. Attempting to hermitize...")
                        H = HamiltonianCpp.hermitize_matrix!(H)  # Use the C++ faithful hermitization
                        if !ishermitian(H)
                            error("Failed to hermitize Hamiltonian matrix. It remains non-Hermitian.")
                        end
                    end
                end
                
                # Track Hamiltonian memory
                H_memory_size = Profiling.get_memory_size(H)
                Profiling.track_memory("Hamiltonian_Matrix", H_memory_size)
                if current_params.profile
                    println("Hamiltonian matrix memory: $(Profiling.format_bytes(H_memory_size))")
                    
                    # Check current memory after Hamiltonian construction
                    current_mem = Base.gc_live_bytes()
                    println("Memory after Hamiltonian construction: $(Profiling.format_bytes(current_mem))")
                end
                
                # Save Hamiltonian matrix with performance tracking
                if args["save-hamiltonian"]
                    Profiling.@time_operation "Hamiltonian_Saving" begin
                        println("Saving Hamiltonian matrix for comparison...")
                        hamiltonian_suffix = hamiltonian_method == "cpp-faithful" ? "hamiltonian_matrix_cpp" : "hamiltonian_matrix_julia"
                        hamiltonian_filename = Utility.create_filename(current_params, 0, hamiltonian_suffix)
                        
                        # Save Hamiltonian in the same format as C++ (row-major, complex numbers)
                        open(hamiltonian_filename, "w") do io
                            println(io, "// $(hamiltonian_method)-generated Hamiltonian matrix")
                            println(io, "// Matrix size: $(size(H, 1))×$(size(H, 2))")
                            println(io, "// Format: row-major order, tab-separated complex numbers per row")
                            println(io, "// Complex format: real+imag*im")
                            
                            for i in 1:size(H, 1)
                                row_entries = String[]
                                for j in 1:size(H, 2)
                                    val = H[i, j]
                                    push!(row_entries, "$(real(val))+$(imag(val))*im")
                                end
                                println(io, join(row_entries, "\t"))
                            end
                        end
                        println("Hamiltonian matrix saved to: $hamiltonian_filename")

                        # Save also the channels for reference
                        channels_suffix = basis_method == "cpp-faithful" ? "channels_cpp" : "channels_julia"
                        channels_filename = Utility.create_filename(current_params, 0, channels_suffix)
                        open(channels_filename, "w") do io
                            println(io, "// $(basis_method)-generated channels")
                            println(io, "// Number of channels: $QN")
                            println(io, "// Format: [j, m, S_actual, ms_actual] for each channel")
                            for ch in channels  
                                println(io, "[", join(ch, ", "), "]")
                            end
                        end
                    end
                end
                
                # Solve the eigenproblem with performance tracking
                println("Starting eigenproblem solution...")
                
                spectra_unit = args["spectra-unit"]
                matrix_size = size(H, 1)
                
                # Prepare eigenproblem solver arguments based on spectra unit
                eig_kwargs = Dict{Symbol, Any}()
                if spectra_unit == "Brot"
                    eig_kwargs[:spectra_in_Brot] = true
                    eig_kwargs[:B_rot_quantity] = current_params.Brot
                elseif spectra_unit == "Hz"
                    eig_kwargs[:spectra_in_Hz] = true
                else # "Hartree" or any other value defaults to atomic units
                    # No conversion needed for Hartree units (default)
                end
                
                eig = Profiling.@time_operation "Eigenvalue_Computation" begin
                    println("INFO: Eigenvalues will be output in $spectra_unit units")
                    Eigenproblem.solve_eigenproblem(H; eig_kwargs...)
                end
                
                solver_suffix = eigensolver_method == "mkl" && MKL_AVAILABLE[] ? "mkl" : "julia"
                
                println("Eigenproblem solved. Found $(length(eig.values)) eigenvalues.")
                
                # Track eigenvalue/eigenvector memory
                eig_memory_size = Profiling.get_memory_size(eig)
                Profiling.track_memory("Eigenvalue_Results", eig_memory_size)
                if current_params.profile
                    println("Eigenvalue results memory: $(Profiling.format_bytes(eig_memory_size))")
                end
                
                # Save results with performance tracking
                Profiling.@time_operation "Results_Saving" begin
                    eigenvalues_suffix = "eigenvalues_$(hamiltonian_method)_$(solver_suffix)_$(spectra_unit)"
                    filename = Utility.create_filename(current_params, 0, eigenvalues_suffix)
                    open(filename, "w") do io
                        # Add header comment indicating units
                        println(io, "# Eigenvalues in $spectra_unit units")
                        println(io, "# Solver: $(eigensolver_method) ($(solver_suffix))")
                        println(io, "# Hamiltonian method: $hamiltonian_method")
                        println(io, "# Matrix size: $(matrix_size)×$(matrix_size)")
                        for val in eig.values
                            println(io, val)
                        end
                    end
                    println("Spectrum calculation complete. Eigenvalues saved to $filename.")
                end
            end
        end
    
    elseif mode == 1
        if vsWhat == 0
            println("Running quench dynamics...")
            flush(stdout)
            Quench.run_quench_dynamics(args)
        elseif vsWhat == 1
            println("Consolidating the quench data from distributed calculations...")
            flush(stdout)
            Quench.consolidate_quench_data(
                current_params; remove_distributed=false
                )
        elseif vsWhat == 2
            println("Plotting the quench dynamics results...")
            flush(stdout)
            Quench.plot_quench_dynamics(
                current_params
                )
        else
            error("Unknown quench option: $vsWhat")
        end
    else
        error("Unknown mode: $mode")
    end
    
    # Print profiling reports at the end
    Profiling.print_profiling_reports()
    Profiling.cleanup_profiling()
end

# This part allows running `julia src/FewMolecules.jl ...`
# It will execute FewMolecules.main()
if abspath(PROGRAM_FILE) == @__FILE__
    # Ensure __init__() is called when the file is run directly
    __init__()
    main()
end

end # module FewMolecules