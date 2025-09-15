module SavingUtils

export save_eigenvalues, save_eigenproblem, save_simulation_params, save_basis

using ..Utility
using ..Parameters
using ..Profiling
using Unitful: ustrip
using Unitful
using LinearAlgebra
using SparseArrays
using Arrow
using JLD2


# Saving functions

# Accept eigenvalue type T instead of hardcoding Float64
function save_eigenvalues(
    eigenvalues::AbstractVector,
    params::Parameters.Params,
    iterator::Int,
    value,
    how_many::Int;
    spectra_unit::String = "Brot",
    subfolder::String = "",
    file_type=".arrow",
)
    Profiling.@time_operation "Saving_Eigenvalues" begin
        filename = Utility.get_filename(params; iterator=iterator, object_name="eigenvalues", subfolder=subfolder, file_extension=file_type)

        if file_type == ".arrow"
            # Build Arrow table with metadata as extra columns
            ham_var_type = Utility.get_var_dtype(params.hamiltonian_dtype)
            tbl = (; (Symbol("Eigenvalues")=>eigenvalues[1:min(how_many, length(eigenvalues))]),)
            Arrow.write(filename, tbl)
            Profiling.track_memory("Saved_Eigenvalues", Profiling.get_memory_size(tbl))
            if !params.progressbar
                if params.print_colored
                    printstyled("    Eigenvalues saved in ", color=:white)
                    printstyled("$(spectra_unit) ", color=:light_blue)
                    printstyled("units to: ", color=:white)
                    printstyled("$filename", color=:magenta)
                    printstyled(" (Arrow format)\n", color=:white)
                else
                    println("    Eigenvalues saved to: $filename (Arrow format)")
                end
            end
        elseif file_type == ".txt"
            # Save as plain text file
            # Save as .txt (human-readable)
            open(filename, "w") do io
                println(io, "# Eigenvalues in $spectra_unit units")
                println(io, "# Iteration value: $value")
                println(io, "# Saved eigenvalues: $(min(how_many, length(eigenvalues)))")
                for i in 1:min(how_many, length(eigenvalues))
                    println(io, eigenvalues[i])
                end
            end
            Profiling.track_memory("Saved_Eigenvalues", Profiling.get_memory_size(eigenvalues))
            if params.print_colored
                printstyled("    Eigenvalues saved to: ", color=:white)
                printstyled("$filename", color=:magenta)
                printstyled(" (txt format)\n", color=:white)
            else
                println("    Eigenvalues saved to: $filename (txt format)")
            end
        else
            error("Unsupported file type: $file_type. Supported types are '.arrow' and '.txt'")
        end
    end
end

function save_eigenproblem(
    eig_result,
    channels::Vector{Vector{Float64}},
    QN_max::Int,
    params::Parameters.Params,
    iterator::Int;
    subfolder::String = "",
    threshold::Float64 = 0.005,
    file_type=".arrow",
)
    Profiling.@time_operation "Saving_Eigenproblem" begin
        values = nothing
        vectors = nothing
        if typeof(eig_result) <: LinearAlgebra.Eigen
            values = eig_result.values
            vectors = eig_result.vectors
        elseif hasfield(typeof(eig_result), :values) && hasfield(typeof(eig_result), :vectors)
            values = eig_result.values
            vectors = eig_result.vectors
        else
            error("save_eigenproblem: eig_result must be LinearAlgebra.Eigen type or as fields .values and .vectors")
        end
        
        filename = Utility.get_filename(params; iterator=iterator, object_name="eigenstates", subfolder=subfolder, file_extension=file_type)

        if file_type == ".arrow"
            # Sparse encoding: each eigenstate is a column of NamedTuples (basis_state_index, coeff_squared)
            n_states = length(values)
            eigenstate_cols = Vector{Vector{NamedTuple{(:basis_state_index, :coeff_squared), Tuple{Int, Float64}}}}(undef, n_states)
            max_len = 0
            empty_tuple = (basis_state_index=-1, coeff_squared=0.0)
            for state in 1:n_states
                col_data = NamedTuple{(:basis_state_index, :coeff_squared), Tuple{Int, Float64}}[]
                for i in 1:QN_max
                    coeff_squared = abs2(vectors[i, state])
                    if coeff_squared > threshold
                        push!(col_data, (basis_state_index=i, coeff_squared=coeff_squared))
                    end
                end
                eigenstate_cols[state] = col_data
                max_len = max(max_len, length(col_data))
            end

            # Pad each column to max_len with missing
            for state in 1:n_states
                while length(eigenstate_cols[state]) < max_len
                    push!(eigenstate_cols[state], empty_tuple)
                end
            end

            # Build Arrow table
            tbl = NamedTuple()
            for state in 1:n_states
                tbl = merge(tbl, (Symbol("eigenstate_$state") => eigenstate_cols[state],))
            end

            Arrow.write(filename, tbl)
            Profiling.track_memory("Saved_Eigenproblem_Arrow", Profiling.get_memory_size(tbl))
            if !params.progressbar
                if params.print_colored
                    printstyled("    Eigenstates saved to: ", color=:white)
                    printstyled("$filename", color=:magenta)
                    printstyled(" (Arrow format)\n", color=:white)
                else
                    println("    Eigenstates saved to: $filename (Arrow format)")
                end
            end
            return dirname(filename)
        elseif file_type == ".txt"
            # Sort eigenvalues and get sorting indices
            sorted_indices = sortperm(values)
            open(filename, "w") do io
                println(io, "# Eigenvectors (eigenstates) for this iteration")
                println(io, "# Columns: coeff^2, followed by quantum numbers for each channel")
                num_states_to_print = length(values)  # Print all states, but can be adjusted if needed
                for state in 1:num_states_to_print
                    idx = sorted_indices[state]
                    energy = values[idx]
                    println(io, "State $state, energy $energy:")
                    for i in 1:QN_max
                        coefficient = vectors[i, idx]
                        coeff_squared = abs2(coefficient)
                        if coeff_squared > 0.005
                            print(io, coeff_squared, "\t")
                            for j in 1:length(channels[i])
                                print(io, channels[i][j], "\t")
                            end
                            println(io)
                        end
                    end
                end
            end
            Profiling.track_memory("Saved_Eigenproblem_Txt", Profiling.get_memory_size(vectors))
            if params.print_colored
                printstyled("    Eigenstates saved to: ", color=:white)
                printstyled("$filename", color=:magenta)
                printstyled(" (txt format)\n", color=:white)
            else
                println("    Eigenstates saved to: $filename (txt format)")
            end
            return dirname(filename)
        elseif file_type == ".jld2"
            vectors[abs2.(vectors) .< threshold] .= 0
            vectors = SparseMatrixCSC(vectors)
            dropzeros!(vectors)
            jldsave(filename; vectors)
            Profiling.track_memory("Saved_Eigenproblem_JLD2", Profiling.get_memory_size(vectors))
            if !params.progressbar
                if params.print_colored
                    printstyled("    Eigenstates saved to: ", color=:white)
                    printstyled("$filename", color=:magenta)
                    printstyled(" (sparse format, .jld2 file)\n", color=:white)
                else
                    println("    Eigenstates saved to: $filename (sparse format, .jld2 file)")
                end
            end
            return dirname(filename)
        else 
            error("Unsupported file type: $file_type. Supported types are '.arrow' and '.txt'")
        end
    end
end

function save_simulation_params(
    params::Parameters.Params;
    )
    Profiling.@time_operation "Saving_Simulation_Params" begin

        filename = "simulation_params_S$(params.spin)_Jmax$(params.JMAX)_Mtot$(params.Mtot).jld2"
        dirname = Utility.get_filename(params; object_name="simulation_params", just_dir=true, override_params=true)
        # Ensure the directory exists
        mkpath(dirname)

        # Create a JLD2 file with the parameters
        jld2_file = joinpath(dirname, filename)
        if !isfile(jld2_file)
            jldsave(jld2_file; params)
            Profiling.track_memory("Saved_Simulation_Params", Profiling.get_memory_size(params))
            if params.print_colored
                printstyled("    Simulation parameters saved to: ", color=:white)
                printstyled("$jld2_file\n", color=:magenta)
            else
                println("    Simulation parameters saved to: $jld2_file")
            end
        else
            if params.print_colored
                printstyled("    Simulation parameters already exist at: ", color=:white)
                printstyled("$jld2_file\n", color=:magenta)
            else
                println("    Simulation parameters already exist at: $jld2_file")
            end
        end
    end
end

function save_basis(
    channels::Vector{Vector{Float64}},
    params::Parameters.Params;
)
    Profiling.@time_operation "Saving_Basis" begin
        filename = "channels_S$(params.spin)_Jmax$(params.JMAX)_Mtot$(params.Mtot).arrow"
        dirname = Utility.get_filename(params; object_name="channels", just_dir=true, override_params=true)
        mkpath(dirname)

        arrow_file = joinpath(dirname, filename)
        if !isfile(arrow_file)
            channels_tbl = (; (Symbol("channel_$i") => channels[i] for i in 1:length(channels))...)
            Arrow.write(arrow_file, channels_tbl)
            Profiling.track_memory("Saved_Basis", Profiling.get_memory_size(channels_tbl))
        end
    end
end

end # module SavingUtils