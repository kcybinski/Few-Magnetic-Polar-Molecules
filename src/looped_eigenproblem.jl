module LoopedEigenproblem

using LinearAlgebra
using SparseArrays
using Unitful
using JSON
using Dates
using ProgressBars
using ..Parameters
using ..Units
using ..Utility
using ..SavingUtils
using ..Eigenproblem
using ..HamiltonianCpp
using ..HamFunctionsCpp
using ..Profiling

export looped_eigenproblem

function looped_eigenproblem(channels::Vector{Vector{Float64}}, QN_max::Int, params::Parameters.Params, how_many_eigenstates_to_save::Int; 
                            param_strings,
                           spectra_unit::String="Brot", 
                           eigensolver_method::String="auto", 
                           mkl_available::Bool=false)
    basis_size = length(channels)
    geometry_length = 2 * (params.how_many_molecules - 1)
    
    println("Starting looped eigenproblem calculation...")
    println("Iterating from $(params.START_IT) to $(params.END_IT) with step $(params.step)")
    println("vsWhat = $(params.vsWhat)")
    println("Eigenvalues will be output in $spectra_unit units")
    
    # Determine which term to skip in base Hamiltonian construction
    skip_terms = Set(Parameters.get_skip_term_for_vsWhat(vsWhat=params.vsWhat))
    if params.skip_terms !== nothing
        skip_terms = union(skip_terms, params.skip_terms)
    end
    println("Building base Hamiltonian, skipping: $skip_terms")

    # Build base Hamiltonian once, excluding the term being iterated
    base_hamiltonian = HamiltonianCpp.buildHamiltonian(
        channels, QN_max, params, true; 
        skip_terms=skip_terms
    )
    
    println("Base Hamiltonian constructed. Starting parameter sweep...")
    
    # Create solver settings to pass to iteration functions
    solver_settings = (
        spectra_unit=spectra_unit,
        eigensolver_method=eigensolver_method,
        mkl_available=mkl_available
    )

    flush(stdout)  # Ensure all output is printed before starting iterations
    
    # Dispatch to appropriate iteration function based on vsWhat
    if params.vsWhat == 0  # Electric field
        println("Iterating over electric field...")
        iterate_electric_field(base_hamiltonian, channels, QN_max, params, param_strings, how_many_eigenstates_to_save, solver_settings)
    elseif params.vsWhat == 1  # Magnetic field
        println("Iterating over magnetic field...")
        iterate_magnetic_field(base_hamiltonian, channels, QN_max, params, param_strings, how_many_eigenstates_to_save, solver_settings)
    elseif params.vsWhat == 2  # Spin-rotation coupling
        println("Iterating over spin-rotation coupling...")
        iterate_spin_rotation(base_hamiltonian, channels, QN_max, params, param_strings, how_many_eigenstates_to_save, solver_settings)
    elseif params.vsWhat == 3  # Geometry/distances
        println("Iterating over geometry scaling...")
        iterate_geometry_scaling(base_hamiltonian, channels, QN_max, params, param_strings, how_many_eigenstates_to_save, solver_settings)
    elseif params.vsWhat == 4  # Electric dipole moment
        println("Iterating over electric dipole moment...")
        iterate_electric_dipole(base_hamiltonian, channels, QN_max, params, param_strings, how_many_eigenstates_to_save, solver_settings)
    elseif params.vsWhat == 5  # Magnetic dipole moment
        println("Iterating over magnetic dipole moment...")
        iterate_magnetic_dipole(base_hamiltonian, channels, QN_max, params, param_strings, how_many_eigenstates_to_save, solver_settings)
    else
        error("Unknown vsWhat value: $(params.vsWhat)")
    end
end

function iterate_electric_field(base_hamiltonian::SparseMatrixCSC, channels::Vector{Vector{Float64}}, QN_max::Int, params::Parameters.Params, param_strings::Vector{String}, how_many_eigenstates_to_save::Int, solver_settings)
    if params.IT_NO_OVERRIDE === nothing
        if params.progressbar
            it_generator = ProgressBar(0:(params.END_IT - params.START_IT))
        else
            it_generator = 0:(params.END_IT - params.START_IT)
        end
    else
        it_generator = params.IT_NO_OVERRIDE
    end

    for iterator in it_generator
        E_value = Utility.get_loop_value(params, iterator)
        if !params.progressbar
            if params.print_colored
                printstyled("Iteration $iterator:", color=:cyan)
                printstyled(" E = $E_value kV/cm", color=:green)
                printstyled(" current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n", color=:white)
            else
                println("Iteration $iterator: E = $E_value kV/cm, current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
            end
        end

        flush(stdout)

        # Create new parameters with updated E field
        new_params = deepcopy(params)
        new_params.E = E_value * 1u"kV/cm"
        Profiling.track_memory("Params_DeepCopy", Profiling.get_memory_size(new_params))
        
        # Build only the electric field term matrix
        electric_field_matrix = HamiltonianCpp.createTermMatrix(channels, QN_max, new_params, :electric_field)
        Profiling.track_memory("Electric_Field_Matrix", Profiling.get_memory_size(electric_field_matrix))
        
        flush(stdout)

        find_eigenproblem_and_save(E_value, electric_field_matrix, iterator, base_hamiltonian, how_many_eigenstates_to_save, new_params, channels, QN_max, solver_settings)
        if params.progressbar
            set_multiline_postfix(it_generator, "E = $E_value kV/cm\nLast iteration at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        end
    end
end

function iterate_magnetic_field(base_hamiltonian::SparseMatrixCSC, channels::Vector{Vector{Float64}}, QN_max::Int, params::Parameters.Params, param_strings::Vector{String}, how_many_eigenstates_to_save::Int, solver_settings)
    if params.IT_NO_OVERRIDE === nothing
        if params.progressbar
            it_generator = ProgressBar(0:(params.END_IT - params.START_IT))
        else
            it_generator = 0:(params.END_IT - params.START_IT)
        end
    else
        start_it_override = params.IT_NO_OVERRIDE - 1
        if params.ITS_PER_OVERRIDE === nothing
            it_generator = start_it_override
        else
            it_generator = start_it_override:(start_it_override + params.ITS_PER_OVERRIDE - 1)
        end
    end
    for iterator in it_generator
        B_value = Utility.get_loop_value(params, iterator)
        if !params.progressbar
            if params.print_colored
                printstyled("Iteration $iterator:", color=:cyan)
                printstyled(" B = $B_value G", color=:green)
                printstyled(" current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n", color=:white)
            else
                println("Iteration $iterator: B = $B_value G, current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
            end
        end

        flush(stdout)

        # Create new parameters with updated B field
        new_params = deepcopy(params)
        new_params.B = B_value * 1u"Gauss"
        Profiling.track_memory("Params_DeepCopy", Profiling.get_memory_size(new_params))

        # Build only the magnetic field term matrix
        magnetic_field_matrix = HamiltonianCpp.createTermMatrix(channels, QN_max, new_params, :magnetic_field)
        Profiling.track_memory("Magnetic_Field_Matrix", Profiling.get_memory_size(magnetic_field_matrix))

        flush(stdout)

        find_eigenproblem_and_save(B_value, magnetic_field_matrix, iterator, base_hamiltonian, how_many_eigenstates_to_save, new_params, channels, QN_max, solver_settings)
        if params.progressbar
            set_multiline_postfix(it_generator, "B = $B_value G\nLast iteration at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        end
    end
end

function iterate_spin_rotation(base_hamiltonian::SparseMatrixCSC, channels::Vector{Vector{Float64}}, QN_max::Int, params::Parameters.Params, param_strings::Vector{String}, how_many_eigenstates_to_save::Int, solver_settings)
    if params.IT_NO_OVERRIDE === nothing
        if params.progressbar
            it_generator = ProgressBar(0:(params.END_IT - params.START_IT))
        else
            it_generator = 0:(params.END_IT - params.START_IT)
        end
    else
        start_it_override = params.IT_NO_OVERRIDE - 1
        if params.ITS_PER_OVERRIDE === nothing
            it_generator = start_it_override
        else
            it_generator = start_it_override:(start_it_override + params.ITS_PER_OVERRIDE - 1)
        end
    end
    for iterator in it_generator
        A_value = Utility.get_loop_value(params, iterator)
        if !params.progressbar
            if params.print_colored
                printstyled("Iteration $iterator:", color=:cyan)
                printstyled(" A = $A_value kHz", color=:green)
                printstyled(" current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n", color=:white)
            else
                println("Iteration $iterator: A = $A_value kHz, current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
            end
        end

        flush(stdout)

        # Create new parameters with updated A
        new_params = deepcopy(params)
        new_params.A = A_value * 1u"kHz"
        Profiling.track_memory("Params_DeepCopy", Profiling.get_memory_size(new_params))

        # Build only the spin-rotation term matrix
        spin_rotation_matrix = HamiltonianCpp.createTermMatrix(channels, QN_max, new_params, :spin_rotation)
        Profiling.track_memory("Spin_Rotation_Matrix", Profiling.get_memory_size(spin_rotation_matrix))

        flush(stdout)

        find_eigenproblem_and_save(A_value, spin_rotation_matrix, iterator, base_hamiltonian, how_many_eigenstates_to_save, new_params, channels, QN_max, solver_settings)
        if params.progressbar
            set_multiline_postfix(it_generator, "A = $A_value kHz\nLast iteration at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        end
    end
end

function iterate_geometry_scaling(base_hamiltonian::SparseMatrixCSC, channels::Vector{Vector{Float64}}, QN_max::Int, params::Parameters.Params, param_strings::Vector{String}, how_many_eigenstates_to_save::Int, solver_settings)
    if params.IT_NO_OVERRIDE === nothing
        if params.progressbar
            it_generator = ProgressBar(0:(params.END_IT - params.START_IT))
        else
            it_generator = 0:(params.END_IT - params.START_IT)
        end
    else
        start_it_override = params.IT_NO_OVERRIDE - 1
        if params.ITS_PER_OVERRIDE === nothing
            it_generator = start_it_override
        else
            it_generator = start_it_override:(start_it_override + params.ITS_PER_OVERRIDE - 1)
        end
    end
    for iterator in it_generator
        R_multiplier = Utility.get_loop_value(params, iterator)
        if !params.progressbar
            if params.print_colored
                printstyled("Iteration $iterator:", color=:cyan)
                printstyled(" R_multiplier = $R_multiplier", color=:green)
                printstyled(" current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n", color=:white)
            else
                println("Iteration $iterator: R_multiplier = $R_multiplier, current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
            end
        end

        flush(stdout)

        # Create new parameters with scaled geometry
        new_params = deepcopy(params)
        # Scale all interaction geometries
        for i in 1:length(new_params.interaction_geometries)
            R_original, angle = new_params.interaction_geometries[i]
            R_scaled = R_original * R_multiplier
            new_params.interaction_geometries[i] = (R_scaled, angle)
        end
        Profiling.track_memory("Params_DeepCopy", Profiling.get_memory_size(new_params))

        # Build dipolar interaction terms with new geometry
        combined_matrix = HamiltonianCpp.createCombinedMatrix(channels, QN_max, new_params, [:electric_dipolar, :magnetic_dipolar])
        Profiling.track_memory("Combined_Interaction_Matrix", Profiling.get_memory_size(combined_matrix))

        flush(stdout)

        find_eigenproblem_and_save(R_multiplier, combined_matrix, iterator, base_hamiltonian, how_many_eigenstates_to_save, new_params, channels, QN_max, solver_settings)
        if params.progressbar
            set_multiline_postfix(it_generator, "R_multiplier = $R_multiplier\nLast iteration at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        end
    end
end

function iterate_electric_dipole(base_hamiltonian::SparseMatrixCSC, channels::Vector{Vector{Float64}}, QN_max::Int, params::Parameters.Params, param_strings::Vector{String}, how_many_eigenstates_to_save::Int, solver_settings)
    if params.IT_NO_OVERRIDE === nothing
        if params.progressbar
            it_generator = ProgressBar(0:(params.END_IT - params.START_IT))
        else
            it_generator = 0:(params.END_IT - params.START_IT)
        end
    else
        start_it_override = params.IT_NO_OVERRIDE - 1
        if params.ITS_PER_OVERRIDE === nothing
            it_generator = start_it_override
        else
            it_generator = start_it_override:(start_it_override + params.ITS_PER_OVERRIDE - 1)
        end
    end
    for iterator in it_generator
        d_el_value = Utility.get_loop_value(params, iterator)
        if !params.progressbar
            if params.print_colored
                printstyled("Iteration $iterator:", color=:cyan)
                printstyled(" d_el = $d_el_value D", color=:green)
                printstyled(" current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n", color=:white)
            else
                println("Iteration $iterator: d_el = $d_el_value D, current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
            end
        end

        flush(stdout)

        # Create new parameters with updated d_el
        new_params = deepcopy(params)
        new_params.d_el = d_el_value * 1u"D"
        Profiling.track_memory("Params_DeepCopy", Profiling.get_memory_size(new_params))

        # Build electric dipolar interaction + electric field terms (as in C++ case 4)
        combined_matrix = HamiltonianCpp.createCombinedMatrix(channels, QN_max, new_params, [:electric_dipolar, :electric_field])
        Profiling.track_memory("Combined_Electric_Dipole_Matrix", Profiling.get_memory_size(combined_matrix))

        flush(stdout)

        find_eigenproblem_and_save(d_el_value, combined_matrix, iterator, base_hamiltonian, how_many_eigenstates_to_save, new_params, channels, QN_max, solver_settings)
        if params.progressbar
            set_multiline_postfix(it_generator, "d_el = $d_el_value D\nLast iteration at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        end
    end
end

function iterate_magnetic_dipole(base_hamiltonian::SparseMatrixCSC, channels::Vector{Vector{Float64}}, QN_max::Int, params::Parameters.Params, param_strings::Vector{String}, how_many_eigenstates_to_save::Int, solver_settings)
    if params.IT_NO_OVERRIDE === nothing
        if params.progressbar
            it_generator = ProgressBar(0:(params.END_IT - params.START_IT))
        else
            it_generator = 0:(params.END_IT - params.START_IT)
        end
    else
        start_it_override = params.IT_NO_OVERRIDE - 1
        if params.ITS_PER_OVERRIDE === nothing
            it_generator = start_it_override
        else
            it_generator = start_it_override:(start_it_override + params.ITS_PER_OVERRIDE - 1)
        end
    end
    for iterator in it_generator
        d_mg_value = Utility.get_loop_value(params, iterator)
        if !params.progressbar
            if params.print_colored
                printstyled("Iteration $iterator:", color=:cyan)
                printstyled(" d_mg = $d_mg_value μB", color=:green)
                printstyled(" current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n", color=:white)
            else
                println("Iteration $iterator: d_mg = $d_mg_value μB, current time = $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
            end
        end

        flush(stdout)

        # Create new parameters with updated d_mg
        new_params = deepcopy(params)
        new_params.d_mg = d_mg_value * 1u"μB"
        Profiling.track_memory("Params_DeepCopy", Profiling.get_memory_size(new_params))

        # Build magnetic dipolar interaction + magnetic field terms (as in C++ case 5)
        combined_matrix = HamiltonianCpp.createCombinedMatrix(channels, QN_max, new_params, [:magnetic_dipolar, :magnetic_field])
        Profiling.track_memory("Combined_Magnetic_Dipole_Matrix", Profiling.get_memory_size(combined_matrix))

        flush(stdout)

        find_eigenproblem_and_save(d_mg_value, combined_matrix, iterator, base_hamiltonian, how_many_eigenstates_to_save, new_params, channels, QN_max, solver_settings)
        if params.progressbar
            set_multiline_postfix(it_generator, "d_mg = $d_mg_value μB\nLast iteration at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        end
    end
end

"""
    find_eigenproblem_and_save(value, matrix_to_add, iterator, param_strings, base_hamiltonian, how_many_eigenstates_to_save, params, channels, QN_max, solver_settings, print_results=false)

Like `find_eigenproblem_and_save`, but always returns the full eigenproblem (eigenvalues and eigenvectors).
"""
function find_eigenproblem_and_save(
    value::Float64, matrix_to_add::SparseMatrixCSC, iterator::Int,
    base_hamiltonian::SparseMatrixCSC, how_many_eigenstates_to_save::Int, params::Parameters.Params,
    channels::Vector{Vector{Float64}}, QN_max::Int, solver_settings, print_results::Bool=false
)
    if !params.progressbar
        if params.print_colored
            printstyled("  Computing eigenproblem for iteration $iterator...", color=:cyan)
        else
            println("  Computing eigenproblem for iteration $iterator...")
        end
    end

    flush(stdout)  # Ensure all output is printed before starting computation

    # Build final Hamiltonian: base + varying term (sparse addition is efficient)
    hamiltonian_sparse = Profiling.@time_operation "Hamiltonian_Construction_Looped" HamFunctionsCpp.hermitize_matrix!(base_hamiltonian + matrix_to_add)
    hamiltonian_dense = Matrix(hamiltonian_sparse)

    # Track memory for Hamiltonian matrix
    Profiling.track_memory("Hamiltonian_Matrix_Looped", Profiling.get_memory_size(hamiltonian_dense))

    # Extract solver settings
    spectra_unit = solver_settings.spectra_unit

    # Prepare eigenproblem solver arguments based on spectra unit
    eig_kwargs = Dict{Symbol, Any}()
    if spectra_unit == "Brot"
        eig_kwargs[:spectra_in_Brot] = true
        eig_kwargs[:B_rot_quantity] = params.Brot
    elseif spectra_unit == "Hz"
        eig_kwargs[:spectra_in_Hz] = true
    end

    if params.mode == 0
        eig_kwargs[:compute_eigenvectors] = false
    elseif params.mode == 2
        eig_kwargs[:compute_eigenvectors] = true
    else
        error("Wrong mode: $(mode)")
    end

    # Always solve full eigenproblem (eigenvalues + eigenvectors)
    eig = Profiling.@time_operation "Eigenvalue_Computation_Looped" begin
        matrix_size = size(hamiltonian_dense, 1)
        if !params.progressbar
            println("    Solving eigenproblem for matrix size: $(matrix_size)×$(matrix_size)")
        end
        Eigenproblem.solve_eigenproblem(hamiltonian_dense; eig_kwargs...)
    end

    Profiling.track_memory("Eigenproblem_Values", Profiling.get_memory_size(eig.values))
    if hasproperty(eig, :vectors)
        Profiling.track_memory("Eigenproblem_Vectors", Profiling.get_memory_size(eig.vectors))
    end

    # Track memory for eigenproblem result
    Profiling.track_memory("Eigenproblem_Result_Looped", Profiling.get_memory_size(eig))

    if print_results
        println("    First 5 eigenvalues ($(spectra_unit) units): ", eig.values[1:min(5, length(eig.values))])
    end

    flush(stdout)  # Ensure all output is printed before saving results

    # Save eigenvalues in 'energies' subfolder
    Profiling.@time_operation "Results_Saving_Looped_Eigenvalues" begin
        SavingUtils.save_eigenvalues(
            eig.values, params, iterator, value, how_many_eigenstates_to_save;
            spectra_unit=spectra_unit, subfolder="energies",
            file_type=".arrow",
        )
    end

    if params.mode == 2
        # Convert eigenvectors to sparse before saving
        vectors_sparse = SparseMatrixCSC(eig.vectors)
        dropzeros!(vectors_sparse)
        values = copy(eig.values)  # Keep a copy of eigenvalues if needed
        eig = nothing
        GC.gc()                # Force garbage collection before saving

        Profiling.track_memory("Eigenproblem_Vectors_Sparse", Profiling.get_memory_size(vectors_sparse))

        Profiling.@time_operation "Results_Saving_Looped_Eigenproblem" begin
            SavingUtils.save_eigenproblem(
                (; values=values, vectors=vectors_sparse), channels, QN_max, params, iterator; subfolder="eigenstates",
                file_type=".jld2",
            )
        end
    end

    # Additional analysis based on mode
    if params.mode == 3 # Magnetization mode
        expected_M = find_expected_magnetization(eig, channels, QN_max)
        save_magnetization(expected_M, params, iterator, value, how_many_eigenstates_to_save)
    end

    flush(stdout)  # Ensure all output is printed before returning

    return eig
end

function find_expected_magnetization(eig, channels::Vector{Vector{Float64}}, QN_max::Int)
    expected_M = Float64[]
    # Sort eigenvalues and get sorting indices
    sorted_indices = sortperm(real(eig.values))
    for state_idx in sorted_indices
        magnetization = 0.0
        for i in 1:QN_max
            coefficient = eig.vectors[i, state_idx]
            coeff_squared = abs2(coefficient)
            # Sum ms values for all molecules (assuming 4 QNs per molecule: j, m, S, ms)
            total_ms = 0.0
            for mol in 1:length(channels[i])÷4
                ms = channels[i][4 + 4*(mol-1)]  # ms is the 4th QN for each molecule
                total_ms += ms
            end
            magnetization += Float64(coeff_squared) * total_ms
        end
        push!(expected_M, magnetization)
    end
    return expected_M
end

function get_iterated_term_folder(param_strings, value)
    # Try to infer the iterated term from the folder structure or vsWhat
    # If param_strings[1] is e.g. "EF11", "MF10", "A100", "D_EL0.1", "D_MG100"
    # Use a mapping based on the folder name or value
    # If not found, fallback to "LOOPED"
    folder_map = Dict(
        0 => "EF$(value)",    # Electric field
        1 => "MF$(value)",    # Magnetic field
        2 => "A$(value)",     # Spin-rotation
        3 => "R$(value)",     # Geometry scaling
        4 => "D_EL$(value)",  # Electric dipole
        5 => "D_MG$(value)"   # Magnetic dipole
    )
    # Try to infer vsWhat from param_strings[1] if possible, else fallback to EF/MF/A/etc.
    # But best is to use the value of vsWhat from the calling context
    # So, parse the folder name from param_strings[1] if possible
    # But for robustness, just use the value and format as above
    # (Assume param_strings[1] is not empty)
    # If value is integer, format without decimal, else with decimal
    if isa(value, Integer)
        value_str = string(value)
    else
        value_str = string(round(value, sigdigits=6))
    end
    # Try to match the folder name to the iterated term
    # If param_strings[1] starts with "EF", "MF", etc., use that
    # Otherwise, fallback to EF/MF/A/etc.
    # For now, just use the mapping
    # To get vsWhat, we can try to parse from param_strings[1], but better to pass vsWhat explicitly if needed
    # For now, fallback to EF/MF/A/etc.
    # If param_strings[1] starts with "EF", use "EF$(value)", etc.
    # But for robustness, just use the mapping for vsWhat=0 (EF), 1 (MF), etc.
    # If param_strings[1] is empty, fallback to "LOOPED"
    # This function is called from within the correct vsWhat context, so we can use a closure or pass vsWhat if needed
    # For now, just use param_strings[1] as a template and replace the number with value
    # But for clarity, just use the mapping for the most common cases
    # If param_strings[1] is not empty, try to match the prefix
    prefix = ""
    if !isempty(param_strings)
        prefix = param_strings[1][1:min(4, end)]
    end
    if startswith(prefix, "EF")
        return "EF$(value_str)"
    elseif startswith(prefix, "MF")
        return "MF$(value_str)"
    elseif startswith(prefix, "A")
        return "A$(value_str)"
    elseif startswith(prefix, "D_EL")
        return "D_EL$(value_str)"
    elseif startswith(prefix, "D_MG")
        return "D_MG$(value_str)"
    elseif startswith(prefix, "R")
        return "R$(value_str)"
    else
        # Fallback: just use EF/MF/A/etc. for vsWhat=0..5
        # This is not perfect, but covers most cases
        # If you want to be more robust, pass vsWhat explicitly
        return "LOOPED"
    end
end

function save_magnetization(magnetization::AbstractVector, params::Parameters.Params, iterator::Int, value, how_many::Int)
    mol_label = Parameters.PARAMS[].how_many_molecules == 1 ? "1_molecule" : "$(Parameters.PARAMS[].how_many_molecules)_molecules"
    override_top_folder = joinpath(mol_label, Utility.get_vs_folder_name(Parameters.PARAMS[].vsWhat))
    filename = Utility.create_filename(params, iterator, "magnetization"; override_top_folder=override_top_folder)
    
    # Ensure the directory exists before writing
    dir_path = dirname(filename)
    if !isdir(dir_path)
        mkpath(dir_path)
    end
    
    open(filename, "w") do io
        print(io, value, "\t")
        for i in 1:min(how_many, length(magnetization))
            print(io, Float64(magnetization[i]), "\t")
        end
        println(io)
    end
    
    println("    Magnetization saved to: $filename")
end

end # module LoopedEigenproblem