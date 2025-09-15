# Quench dynamics logic - faithful port from C++ implementation
module Quench

export run_quench_dynamics

using LinearAlgebra
using SparseArrays
using Unitful
using Printf
using ProgressBars
using Dates
using Base.Threads

using ..Parameters
using ..Basis
using ..BasisCppFaithful  # Use C++ faithful basis generation
using ..HamiltonianCpp  # Use C++ faithful implementation
using ..Eigenproblem
using ..Utility
using ..Units
using ..SavingUtils

"""
Validate that energy and time units are consistent for harmonic evolution.
"""
function validate_harmonic_units_consistency(harmonic_freq, time_harmonic, energy_harmonic)
    println("Harmonic units consistency check:")
    println("  Harmonic frequency: $(harmonic_freq)")
    println("  Time in harmonic units: $(time_harmonic) × (2π/ω)")
    println("  Energy in harmonic units: $(energy_harmonic) × ℏω")
    println("  Phase evolution: exp(-i × $(energy_harmonic) × $(time_harmonic) × 2π)")
    
    # Check if phase makes physical sense
    typical_phase = abs(energy_harmonic * time_harmonic * 2π)
    if typical_phase > 100
        @warn "Large phase evolution detected ($(typical_phase)). Consider adjusting time step or units."
    elseif typical_phase < 0.01
        @warn "Very small phase evolution detected ($(typical_phase)). Time step might be too small or energy scale too small."
    else
        println("  ✓ Phase evolution scale looks reasonable: $(typical_phase)")
    end
end

"""
Faithful port of C++ chosenInitialStateBasisCoefficients for complex matrices.
Returns coefficients of the chosen eigenstate (default: ground state) in the basis.
"""
function chosen_initial_state_basis_coefficients(H::AbstractMatrix, basis_size::Int, which_state::Int=0; args=Dict())
    # Solve eigenproblem for initial Hamiltonian using centralized auto-selection
    println("Solving initial Hamiltonian eigenproblem for chosen state, Eigenstate $which_state...")
    eig_result = Eigenproblem.solve_eigenproblem(
        H,
    )

    # Get the index of the chosen state (1 = ground state in Julia indexing)
    chosen_index = which_state

    # Extract coefficients for the chosen state
    chosen_coefficients = ComplexF64[]
    for i in 1:basis_size
        push!(chosen_coefficients, eig_result.vectors[i, chosen_index])
    end

    println("Coefficients C_{phi,alpha} of the chosen state of H_0 in the basis alpha are noted.")
    return chosen_coefficients
end

"""
Faithful port of C++ expressProductStateInBasis.
Express a product state in the computational basis.
"""
function express_product_state_in_basis(state::Vector{Float64}, channels::Vector{Vector{Float64}})
    QN_per_mol = 4  # j, m, S, ms per molecule
    how_many_molecules = length(state) ÷ QN_per_mol
    basis_size = length(channels)
    
    if length(state) != how_many_molecules * QN_per_mol
        error("Product state disagrees with the number of molecules in the system!")
    end
    
    basis_coefficients = ComplexF64[]
    any_basis_state_found = false
    
    for i in 1:basis_size
        # Check if this basis state matches the product state
        channel = channels[i]
        condition = true
        
        for qn_idx in 1:length(state)
            if abs(channel[qn_idx] - state[qn_idx]) > 1e-12
                condition = false
                break
            end
        end
        
        if condition
            any_basis_state_found = true
            println("Found matching basis state at index $i: $channel")
            push!(basis_coefficients, ComplexF64(1.0, 0.0))  # Coefficient is 1 for the matching state
        else
            push!(basis_coefficients, ComplexF64(0.0, 0.0))
        end
    end
    
    if any_basis_state_found
        println("Coefficients of the product state in the basis alpha are noted.")
    else
        error("Chosen product state is impossible for the chosen parameters. Check Mtot.")
    end
    
    return basis_coefficients
end

"""
Optimized version with precomputed observables and vectorized operations.
"""
function calculate_phi_psi_j_overlaps(eig_final, chosen_coefficients::Vector{ComplexF64}, 
                                               basis_size::Int, save_overlaps::Bool=false)
    
    # Vectorized computation using matrix-vector multiplication
    # This is much faster than the loop-based approach
    phi_j_overlap = eig_final.vectors' * chosen_coefficients
    
    # Check for non-zero imaginary parts (vectorized)
    imag_parts = imag.(phi_j_overlap)
    large_imag = findall(x -> abs(x) > 1e-12, imag_parts)
    
    if !isempty(large_imag)
        for idx in large_imag
            println("Warning: Non-zero imaginary part in overlap for state $idx: $(imag_parts[idx])")
        end
    end
    
    if save_overlaps
        println("Overlap saving functionality can be implemented if needed")
    end
    
    return phi_j_overlap
end

"""
Calculate matrix elements <Ψⱼ|O|ψⱼ'> for given observable.
Highly optimized matrix element calculation using BLAS operations.
"""
function calculate_matrix_of_psi_o_psi_overlaps(eig_final, channels::Vector{Vector{Float64}}, 
                                                        which_observable::String, 
                                                        basis_size::Int;
                                                        use_sparse::Bool=true, 
                                                        obs_matrix=nothing,
                                                        upper_triangle_only::Bool=true,
                                                        )
    
    # Get or compute observable values for all basis states
    if obs_matrix === nothing
        obs_values = zeros(Float64, basis_size)
        QN_per_mol = 4
        how_many_molecules = length(channels[1]) ÷ QN_per_mol
        
        @inbounds for i in 1:basis_size
            obs_val = 0.0
            
            if which_observable == "Mg" || which_observable == "Ms"
                for mol in 1:how_many_molecules
                    ms_idx = (mol-1) * QN_per_mol + 4
                    obs_val += channels[i][ms_idx]
                end
            elseif which_observable == "M"
                for mol in 1:how_many_molecules
                    m_idx = (mol-1) * QN_per_mol + 2
                    obs_val += channels[i][m_idx]
                end
            elseif startswith(which_observable, "ms")
                mol_num = parse(Int, which_observable[3:end])
                if mol_num <= how_many_molecules
                    ms_idx = (mol_num-1) * QN_per_mol + 4
                    obs_val = channels[i][ms_idx]
                end
            elseif startswith(which_observable, "m")
                mol_num = parse(Int, which_observable[2:end])
                if mol_num <= how_many_molecules
                    m_idx = (mol_num-1) * QN_per_mol + 2
                    obs_val = channels[i][m_idx]
                end
            else
                error("Unknown observable: $which_observable")
            end
            
            obs_values[i] = obs_val
        end
    else
        # Use precomputed values
        obs_idx = findfirst(x -> x == which_observable, ["M", "Ms", "m1", "ms1", "m2", "ms2", "m3", "ms3"])
        obs_values = obs_matrix[:, obs_idx]
    end
    
    # Method optimized for MKL/Accelerate
    # Use in-place operations to minimize memory allocations
    V = eig_final.vectors  # Alias for readability
    
    # Step 1: Scale columns of V by observable values (V_scaled = V * diag(obs))
    V_scaled = similar(V)
    for j in 1:basis_size
        for i in 1:basis_size
            V_scaled[i, j] = V[i, j] * obs_values[i]
        end
    end
    
    # Step 2: Compute V† * V_scaled using optimized GEMM
    # This is the bottleneck operation that benefits most from MKL/Accelerate
    meanvalues_complex = V' * V_scaled 

    # Check for imaginary parts
    imag_check = imag.(meanvalues_complex)
    if maximum(abs.(imag_check)) > 1e-12
        println("Warning: Non-zero imaginary parts detected in observable matrix elements!")

    end
    
    println("Matrix of <Ψ|$(which_observable)|Ψ> overlaps calculated using optimized method.")
    return upper_triangle_only ? triu(real.(meanvalues_complex)) : real.(meanvalues_complex)
end

"""
Precompute all observable values for all basis states to avoid redundant calculations in parallel execution.
"""
function precompute_observables(channels::Vector{Vector{Float64}}, observables::Vector{String})
    QN_per_mol = 4
    how_many_molecules = length(channels[1]) ÷ QN_per_mol
    basis_size = length(channels)
    
    # Preallocate matrix for all observables
    obs_matrix = zeros(Float64, basis_size, length(observables))
    
    for (obs_idx, observable) in enumerate(observables)
        @inbounds for i in 1:basis_size
            obs_val = 0.0
            
            if observable == "Mg" || observable == "Ms"
                # Total magnetic moment - vectorized sum
                for mol in 1:how_many_molecules
                    ms_idx = (mol-1) * QN_per_mol + 4
                    obs_val += channels[i][ms_idx]
                end
            elseif observable == "M"
                # Total orbital angular momentum projection
                for mol in 1:how_many_molecules
                    m_idx = (mol-1) * QN_per_mol + 2
                    obs_val += channels[i][m_idx]
                end
            elseif startswith(observable, "ms")
                # Specific molecule ms
                mol_num = parse(Int, observable[3:end])
                if mol_num <= how_many_molecules
                    ms_idx = (mol_num-1) * QN_per_mol + 4
                    obs_val = channels[i][ms_idx]
                end
            elseif startswith(observable, "m")
                # Specific molecule m
                mol_num = parse(Int, observable[2:end])
                if mol_num <= how_many_molecules
                    m_idx = (mol_num-1) * QN_per_mol + 2
                    obs_val = channels[i][m_idx]
                end
            else
                error("Unknown observable: $observable")
            end
            
            obs_matrix[i, obs_idx] = obs_val
        end
    end
    
    return obs_matrix
end

"""
Parallel version for multiple observables computation.
"""
function calculate_all_observables_parallel(eig_final, channels::Vector{Vector{Float64}}, 
                                          observables::Vector{String}, basis_size::Int;
                                          use_precomputed::Bool=true,
                                          )
    
    results = Dict{String, Matrix{Float64}}()
    
    if use_precomputed
        println("Precomputing all observable values...")
        obs_matrix = precompute_observables(channels, observables)
        println("Precomputation complete.")
        
        # Parallel computation of all observables
        @threads for (idx, observable) in collect(enumerate(observables))
            println("Thread $idx started computing observable $observable...")
            flush(stdout)

            results[observable] = calculate_matrix_of_psi_o_psi_overlaps(
                eig_final, channels, observable, basis_size; obs_matrix,
            )
        end
    else
        # Sequential with individual optimization
        for observable in observables
            println("Started computing observable $observable...")
            results[observable] = calculate_matrix_of_psi_o_psi_overlaps(
                eig_final, channels, observable, basis_size
            )
        end
    end
    
    return results
end

function get_quench_time_iterator(params::Parameters.Params; double_counting::Bool=false)
    tmax_harmonic = params.tmax  # dimensionless
    dt_harmonic = params.dt      # dimensionless
    time_steps = Int(ceil(tmax_harmonic / dt_harmonic))

    if params.progressbar
        time_iterator = ProgressBar(0:time_steps)
    elseif params.IT_NO_OVERRIDE !== nothing
        if params.ITS_PER_OVERRIDE !== nothing
            if !double_counting
                time_iterator = params.IT_NO_OVERRIDE:params.IT_NO_OVERRIDE + params.ITS_PER_OVERRIDE - 1
            else 
                time_iterator = params.IT_NO_OVERRIDE:params.IT_NO_OVERRIDE + params.ITS_PER_OVERRIDE
            end
        else
            time_iterator = params.IT_NO_OVERRIDE # Single iteration for distributed runs
        end
    else
        time_iterator = 0:time_steps
    end

    return time_iterator
end

"""
Main quench dynamics calculation.
Faithful port of C++ quenchDynamics function.
Time is handled in atomic units throughout, consistent with C++ implementation.
"""
function quench_dynamics_calculation(H_final::AbstractMatrix, channels::Vector{Vector{Float64}}, 
                                   chosen_coefficients::Vector{ComplexF64}, 
                                   save_couplings::Bool=false; args::Dict, params::Parameters.Params)
    basis_size = length(channels)

    current_params = Parameters.PARAMS[]

    # Check if harmonic frequency is provided and convert time parameters accordingly
    # Get harmonic frequency and time parameters
    if current_params.harmonic_frequency === nothing
        error("Harmonic frequency not set in parameters")
    end

    harmonic_freq = current_params.harmonic_frequency
    println("Using harmonic frequency: $(harmonic_freq) for time evolution")

    flush(stdout)
    # Solve final Hamiltonian eigenproblem using centralized auto-selection. Energies will be in atomic units (Hartree).
    println("Solving final Hamiltonian eigenproblem...")
    eig_final = Eigenproblem.solve_eigenproblem(
        H_final,
    )

    # Convert eigenvalues to harmonic frequency units for consistency with time
    println("Converting eigenvalues to harmonic frequency units...")
    eig_final_harmonic = merge(eig_final, (; values = Units.eigenvalues_au_to_harmonic_units(eig_final.values, harmonic_freq)))

    println("Energy scale comparison:")
    println("  First eigenvalue (Hartree): $(eig_final.values[1])")
    println("  First eigenvalue (harmonic units): $(eig_final_harmonic.values[1])")
    println("  Energy unit: ℏω where ω = 2π × $(harmonic_freq / 2π)")

    flush(stdout)
    # Calculate overlaps between initial state and final eigenstates
    println("Calculating overlaps...")
    phi_j_overlap = calculate_phi_psi_j_overlaps(eig_final, chosen_coefficients, basis_size, true)


    # Calculate matrix elements for each observable
    meanvalues_dict = calculate_all_observables_parallel(eig_final, channels, current_params.observables, basis_size)
    
    # Time parameters are now in harmonic units (dimensionless)
    tmax_harmonic = current_params.tmax  # dimensionless
    dt_harmonic = current_params.dt      # dimensionless

    # Calculate time steps
    time_steps = Int(ceil(tmax_harmonic / dt_harmonic))

    println("Time evolution in harmonic units:")
    println("  tmax = $(tmax_harmonic) × (2π/ω)")
    println("  dt = $(dt_harmonic) × (2π/ω)")
    println("  energy differences in units of ℏω")
    println("  steps = $time_steps")
    println("  characteristic period = $(1.0) × (2π/ω), in SI units: $(uconvert(u"s", 1.0 / (harmonic_freq/2π)))/2π")

    # Storage for time evolution results
    results_dict = Dict{String, Vector{Float64}}()
    for obs in current_params.observables
        results_dict[obs] = Float64[]
    end

    relevant_couplings = Dict{String, Any}()
    if save_couplings
        for obs in current_params.observables
            relevant_couplings[obs] = Dict(
                "js" => Int[], "jprims" => Int[], "amplitudes" => Float64[]
            )
        end
    end

    flush(stdout)

    # Sanity check for time evolution parameters
    # validate_harmonic_units_consistency(harmonic_freq, dt_harmonic, eig_final_harmonic.values[1] - eig_final_harmonic.values[1])

    time_iterator = get_quench_time_iterator(current_params)

    println("IT_NO_OVERRIDE: $(params.IT_NO_OVERRIDE), ITS_PER_OVERRIDE: $(params.ITS_PER_OVERRIDE)")
    println("  iterations in thread = $(length(time_iterator))")
    # Time evolution loop
    for step in time_iterator
        time_harmonic = step * dt_harmonic  # Time in atomic units (dimensionless)
        
        for obs in current_params.observables
            observable_val = 0.0
            meanvalues = meanvalues_dict[obs]
            
            for j in 1:basis_size
                for jprim in 1:basis_size
                    if j == jprim
                        # Diagonal terms
                        mean_val = meanvalues[j, j]
                        observable_val += abs(phi_j_overlap[j])^2 * mean_val
                        
                    elseif j < jprim
                        # Off-diagonal terms with time evolution
                        mean_val = meanvalues[j, jprim]

                        # Convert energy difference to harmonic units
                        energy_diff_harmonic = eig_final_harmonic.values[j] - eig_final_harmonic.values[jprim]
                        
                        overlaps_product = phi_j_overlap[j] * conj(phi_j_overlap[jprim])

                        # Time evolution: exp(-i * ΔE_harmonic * t_harmonic)
                        # exp(-i * ΔE_harmonic * t_harmonic)
                        # Both ΔE and t are now in consistent harmonic units
                        time_factor = exp(-1im * energy_diff_harmonic * time_harmonic)
                        complex_amplitude = overlaps_product * time_factor
                        
                        amplitude = 2 * real(overlaps_product) * mean_val
                        observable_val += 2 * real(complex_amplitude) * mean_val
                        
                        # Save significant couplings
                        if abs(amplitude) > 0.0005 && time_harmonic == 0.0 && save_couplings
                            push!(relevant_couplings[obs]["js"], j)
                            push!(relevant_couplings[obs]["jprims"], jprim)
                            push!(relevant_couplings[obs]["amplitudes"], abs(amplitude))
                        end
                    end
                end
            end
            
            push!(results_dict[obs], observable_val)
        end
        
        if !params.progressbar
            if params.ITS_PER_OVERRIDE === nothing && params.IT_NO_OVERRIDE === nothing
                print_cond = (step % 100 == 0) ? time_steps > 100 : (step % 5 == 0)
            else
                print_cond = true
            end
            if print_cond
                println("t = $(time_harmonic) × (2π/ω) / $(tmax_harmonic) × (2π/ω)")
                obs_vals = join(["$(obs)=$(@sprintf("%.3e", results_dict[obs][end]))" for obs in current_params.observables], ", ")
                println("t = $(time_harmonic) × (2π/ω): $obs_vals")
                flush(stdout)
            end
        else
            obs_vals = join(["$(obs)=$(@sprintf("%.3e", results_dict[obs][end]))" for obs in current_params.observables], ", ")
            set_multiline_postfix(time_iterator, "t = $(time_harmonic) × (2π/ω) / $(tmax_harmonic) × (2π/ω)\n$obs_vals\nLast iteration at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        end
    end
    
    return results_dict, relevant_couplings
end

"""
Main quench dynamics runner function.
"""
function run_quench_dynamics(args::Dict)
    println("Starting quench dynamics calculation...")
    
    # Get current parameters
    current_params = Parameters.PARAMS[]
    
    # Setup for quench
    basis_method = get(args, "basis-method", "cpp-faithful")
    hamiltonian_method = get(args, "hamiltonian-method", "cpp-faithful")
    
    # Generate basis
    println("Generating basis...")
    
    if basis_method == "cpp-faithful"
        channels = BasisCppFaithful.generate_channels_cpp_style(
            current_params.Mtot, false;
            JMAX=current_params.JMAX,
            spin=current_params.spin / 2.0,
            geometry=args["geometry..."]  # Use CLI geometry directly, same as spectrum
        )
    else
        error("Only cpp-faithful basis method is currently supported for quench dynamics")
    end
    
    basis_size = length(channels)
    println("Basis size: $basis_size")

    flush(stdout)
    
    if basis_size == 0
        error("No channels generated. Check JMAX, spin, Mtot, and geometry parameters.")
    end
    
    
    # Get initial state coefficients
    println("Determining initial state...")
    chosen_coefficients = ComplexF64[]
    which_state = get(args, "initial-eigenstate", 1)  # Default to ground state (1 in Julia)

    productstate = current_params.productstate
    current_params.initial_state_label = ""
    if productstate !== nothing
        # Parse product state from string
        state_values = [parse(Float64, strip(x)) for x in split(productstate, ",")]
        chosen_coefficients = express_product_state_in_basis(state_values, channels)
        println("Using product state evolution with state: $state_values")
        current_params.initial_state_label = "Product State [" * join(state_values, ";") * "]"
    elseif current_params.initial_basis_state !== nothing
        which_state = current_params.initial_basis_state
        chosen_coefficients = zeros(ComplexF64, basis_size)
        chosen_coefficients[which_state] = ComplexF64(1.0, 0.0)  # Set the chosen state to 1
        println("Using product state evolution with basis state No. : $(channels[which_state])")
            current_params.initial_state_label = "Product State [" * join(channels[which_state], ";") * "]"
    else
        # Build initial Hamiltonian (H_0)
        println("Building initial Hamiltonian...")
        if hamiltonian_method == "cpp-faithful"
            H_initial = HamiltonianCpp.buildHamiltonian(channels, basis_size, current_params)
        else
            error("Only cpp-faithful Hamiltonian method is currently supported for quench dynamics")
        end
        # Use a specified eigenstate of initial Hamiltonian
        # The default is the ground state (0)
        H_initial_dense = Matrix(H_initial)
        chosen_coefficients = chosen_initial_state_basis_coefficients(H_initial_dense, basis_size, which_state; args=args)
        println("Using eigenstate $which_state of initial Hamiltonian")
        current_params.initial_state_label = "Eigenstate $which_state"
    end

    SavingUtils.save_basis(channels, current_params)
    SavingUtils.save_simulation_params(current_params)
    
    flush(stdout)
    # Build final Hamiltonian (H_0 + quench fields) by regenerating with modified parameters
    println("Building final Hamiltonian with quench fields...")
    
    # Create modified parameters for quench
    quench_params = deepcopy(current_params)
    
    
    quench_applied = false
    if current_params.quenched_E !== nothing
        println("Applying quenched electric field: $(ustrip(current_params.quenched_E)) kV/cm")
        quench_params.E = ustrip(current_params.quenched_E) * u"kV/cm"
        quench_applied = true
    else
        println("No electric field quench - using initial field: $(ustrip(current_params.E)) kV/cm")
    end
    
    if current_params.quenched_B !== nothing
        println("Applying quenched magnetic field: $(ustrip(current_params.quenched_B)) G")
        quench_params.B = ustrip(current_params.quenched_B) * u"Gauss"
        quench_applied = true
    else
        println("No magnetic field quench - using initial field: $(ustrip(current_params.B)) G")
    end
    
    if current_params.quenched_R !== nothing
        println("Applying quenched distance: $(join([(current_params.quenched_R...)...], ", "))")
        
        quench_params.interaction_geometries = current_params.quenched_R
        quench_applied = true
    else
        println("No distance quench - using initial geometry")
    end
    
    if !quench_applied
        @warn "No quench parameters specified! This will be identical to initial state evolution."
        println("Available quench parameters: quenched_E, quenched_B, quenched_R")
    end

    flush(stdout)

    if !isempty(quench_params.skip_terms)
        println("Applying quench with skip terms: $(quench_params.skip_terms)")
    end

    # Regenerate Hamiltonian with quench parameters
    kwargs = Dict{Symbol, Any}()
    if !isempty(quench_params.skip_terms)
        kwargs[:skip_terms] = quench_params.skip_terms
    end
    H_final = HamiltonianCpp.buildHamiltonian(channels, basis_size, quench_params; kwargs...)

    # Convert to dense matrix for time evolution
    H_final_dense = Matrix(H_final)

    # Run quench dynamics
    println("Running quench dynamics time evolution...")
    flush(stdout)

    # Pass observables list to calculation
    results, couplings = quench_dynamics_calculation(
        H_final_dense, channels, chosen_coefficients,
        true; # save_couplings
        args=args,  # Pass args for unit conversion and other parameters
        params=current_params,
    )

    # Save results
    println("Saving quench dynamics results...")
    flush(stdout)

    time_iterator = get_quench_time_iterator(current_params)

    # Save time evolution data
    for (observable, time_data) in results
        filename = Utility.get_filename(current_params, object_name=observable, ensure_mkdir=true)

        open(filename, "w") do f
            # CSV header with parameter info as comments
            if current_params.IT_NO_OVERRIDE === nothing || current_params.IT_NO_OVERRIDE == 0
                println(f, "# Quench dynamics: $observable vs time")
                println(f, "# Time units: harmonic (2π/ω) where ω = 2π × $(current_params.harmonic_frequency/(2π))")
                println(f, "# Initial state: $(current_params.initial_state_label)")
                println(f, "# tmax = $(current_params.tmax) × (2π/ω)")
                println(f, "# dt = $(current_params.dt) × (2π/ω)")
                println(f, "# Initial conditions: E=$(ustrip(current_params.E)) kV/cm, B=$(ustrip(current_params.B)) G")
                if current_params.quenched_E !== nothing
                    println(f, "# quenched_E=$(ustrip(current_params.quenched_E)) kV/cm")
                end
                if current_params.quenched_B !== nothing
                    println(f, "# quenched_B=$(ustrip(current_params.quenched_B)) G")
                end
                if current_params.quenched_R !== nothing
                    println(f, "# quenched_R=$(join([(current_params.quenched_R...)...], " "))")
                end
                skipping_list = current_params.skip_terms !== nothing ? current_params.skip_terms : String[]
                println(f, "# Skipping terms: ", join(skipping_list, "; ", " and "))
                # CSV column header
                println(f, "time_harmonic,value")
            end
            for (time_tmp, value) in zip(time_iterator, time_data)
                time_harmonic = (time_tmp-1) * current_params.dt # Time in harmonic units
                println(f, "$time_harmonic,$value")
            end
        end
        
        println("Saved $observable evolution to: $filename")
        flush(stdout)
    end

    # Create quench identifier string for filenames
    quench_identifier = String[]
    if current_params.quenched_E !== nothing
        push!(quench_identifier, "Eq$(ustrip(current_params.quenched_E))")
    end
    if current_params.quenched_B !== nothing
        push!(quench_identifier, "Bq$(ustrip(current_params.quenched_B))")
    end
    if current_params.quenched_R !== nothing
        push!(quench_identifier, "Rq$(join(ustrip.([(current_params.quenched_R...)...]), "_"))")
    end
    quench_suffix = isempty(quench_identifier) ? "no_quench" : join(quench_identifier, "_")
    
    if current_params.IT_NO_OVERRIDE === nothing || current_params.IT_NO_OVERRIDE == 0
        # Save coupling information
        for (observable, coupling_data) in couplings
            if !isempty(coupling_data["js"])
                filename = Utility.get_filename(current_params, object_name=observable, ensure_mkdir=true, couplings=true)
                open(filename, "w") do f
                    # CSV header with parameter info as comments
                    println(f, "# Significant couplings for $observable")
                    println(f, "# Initial state: $(current_params.initial_state_label)")
                    println(f, "# Quench: $quench_suffix")
                    skipping_list = current_params.skip_terms !== nothing ? current_params.skip_terms : String[]
                    println("Skipping terms: ", join(skipping_list, "; ", " and "))
                    println(f, "# Format: j,jprim,amplitude")
                    println(f, "j,jprim,amplitude")
                    for i in 1:length(coupling_data["js"])
                        j = coupling_data["js"][i]
                        jprim = coupling_data["jprims"][i]
                        amplitude = coupling_data["amplitudes"][i]
                        println(f, "$j,$jprim,$amplitude")
                    end
                end
                println("Saved $observable couplings to: $filename")
            end
        end
    end
    
    println("Quench dynamics calculation completed successfully!")
    return results
end

"""
This function tidies up the distributed quench dynamics results, merging the 
individual results into a consolidated format, and removing any temporary files.
"""
function consolidate_quench_data(params::Parameters.Params; counter_double_counting = false, overwrite_times = true, remove_distributed = true)
    println("Determining initial state...")
    chosen_coefficients = ComplexF64[]
    which_state = current_params.initial_eigenstate
    
    productstate_str = current_params.productstate
    current_params.initial_state_label = ""
    if productstate_str !== nothing && !isempty(productstate_str)
        # Parse product state from string
        state_values = [parse(Float64, strip(x)) for x in split(productstate_str, ",")]
        println("Using product state evolution with state: $state_values")
        current_params.initial_state_label = "Product State [" * join(state_values, ",") * "]"
    else
        # Use a specified eigenstate of initial Hamiltonian
        println("Using eigenstate $which_state of initial Hamiltonian")
        current_params.initial_state_label = "Eigenstate $which_state"
    end
    
    num_iterations = Int(rationalize(params.tmax / params.dt))
    data_folder = Utility.get_filename(params, override=true, just_dir=true, save_dir=params.save_root_folder)

    # # List all observable files
    list_of_observable_files = Glob.glob(joinpath(data_folder, "*.csv"))

    all_times = [it*params.dt for it in 0:num_iterations-1]

    obs_val_dict = Dict([(obs => fill(NaN, num_iterations)) for obs in params.observables])
    obs_time_dict = Dict([(obs => fill(NaN, num_iterations)) for obs in params.observables])
    comment_dict = Dict([(obs => String[]) for obs in params.observables])

    for filename_tmp in list_of_observable_files
        println("Processing file: $filename_tmp")
        if occursin("couplings", filename_tmp)
            println("   Skipping coupling file: $filename_tmp")
            continue  # Skip coupling files
        end
        pattern_low = r"IT\d+-"
        pattern_high = r"-\d+_"
        pattern_obs = r"_[A-Za-z]+\d*.csv"
        m_low = match(pattern_low, basename(filename_tmp))
        m_high = match(pattern_high, basename(filename_tmp))
        it_low_bd = parse(Int, m_low.match[3:end-1])
        it_high_bd = parse(Int, m_high.match[2:end-1])
        obs = match(pattern_obs, basename(filename_tmp)).match[2:end-4]

        if occursin("IT0", basename(filename_tmp))
            open(filename_tmp, "r") do io
                for line in eachline(io)
                    startswith(line, "#") && push!(comment_dict[obs], line)
                end
            end
            csv_tmp = CSV.read(filename_tmp, DataFrame; comment="#", header=0, skipto=9)
        else
            csv_tmp = CSV.read(filename_tmp, DataFrame; comment="#", header=0)
        end

        if counter_double_counting
            it_high_bd -= 1  # Adjust for double counting
        end

        for (it_fnum, itnum_tmp) in enumerate(it_low_bd:it_high_bd)
            obs_val_dict[obs][itnum_tmp + 1] = csv_tmp.Column2[it_fnum]
            if overwrite_times
                obs_time_dict[obs][itnum_tmp + 1] = all_times[itnum_tmp + 1]
            else
                obs_time_dict[obs][itnum_tmp + 1] = csv_tmp.Column1[it_fnum]
            end
        end
    end

    for (idx, obs) in enumerate(params.observables)
        num_nans = count(isnan, obs_val_dict[obs])
        if idx == 1
            println("Number of errored jobs: $(num_nans/params.ITS_PER_OVERRIDE)")
        end

        println("Number of NaNs in obs_val_dict[\"$obs\"]: $num_nans")
    end

    dfs_observables = Dict([
        (obs => DataFrame()) for obs in params.observables
    ])

    for obs in params.observables
        dfs_observables[obs].time_harmonic = obs_time_dict[obs]
        dfs_observables[obs].values = obs_val_dict[obs]
    end

    tmp_dir = mktempdir(data_folder; prefix="consolidated_quench_data__")

    for obs in params.observables
        csv_file = joinpath(tmp_dir, "quench_observable_$(obs).csv")
        # First write the comments to the CSV file
        open(csv_file, "w") do io
            for comment in comment_dict[obs]
                println(io, comment)
            end
        end
        # Then append the DataFrame to the CSV file
        CSV.write(csv_file, dfs_observables[obs]; header=["time_harmonic", "values"], append=true)
        println("Wrote consolidated data for $obs to $csv_file")
    end

    if remove_distributed
        # Remove all the CSV files with `IT` in name from the parent folder
        # This is to clean up the folder after consolidation
        for filename_tmp in list_of_observable_files
            if !occursin("IT", filename_tmp)
                continue  # Skip files that do not match the pattern
            end
            rm(filename_tmp)
        end
    end
end

end # module