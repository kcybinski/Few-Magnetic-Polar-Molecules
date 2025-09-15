# Eigenproblem/diagonalization logic will go here.
module Eigenproblem

export solve_eigenproblem
using LinearAlgebra
using ..Units # For unit conversions
using ..Parameters # For accessing B_rot if needed
using Unitful

# Dense matrix memory estimates
# https://julialinearalgebra.github.io/BLASBenchmarksCPU.jl/dev/memory-required/

# Simple interface using Julia's built-in eigen solver - now with in-place diagonalization
function solve_eigenproblem(H::AbstractMatrix; 
                           spectra_in_Brot::Bool=false,
                           spectra_in_Hz::Bool=false,
                           B_rot_quantity=nothing,
                           compute_eigenvectors::Bool=true,
                           energies_window=nothing,
                           debug_prints=false
                           )
    
    # Handle unit conversions directly on the input matrix (in-place)
    if spectra_in_Brot && B_rot_quantity !== nothing
        B_rot_au = Units.convert_param_to_au(B_rot_quantity, :Brot)
        H ./= ustrip(B_rot_au)  # In-place division
        if debug_prints
            println("Converting eigenvalues to B_rot units (B_rot = $(ustrip(B_rot_au)) Hartree)")
        end
    elseif spectra_in_Hz
        # Convert from Hartree to Hz
        hartree_to_hz = ustrip(auconvert(u"Hz", 1))
        H .*= hartree_to_hz  # In-place multiplication
        if debug_prints
            println("Converting eigenvalues to Hz units (1 Hartree = $(hartree_to_hz) Hz)")
        end
    else
        if debug_prints
            println("Eigenvalues will be in Hartree (atomic units)")
        end
    end
    
    if debug_prints
        # Use Julia's eigen solver directly on H (which may have been modified in-place)
        println("Starting to solve the eigenproblem using Julia built-in solver.")
        println("Using backend : ", LinearAlgebra.BLAS.get_config())
    end
    
    if compute_eigenvectors
        result = eigen!(H)
        values = result.values
        standardized_vectors = standardize_eigenvector_phases(result.vectors)
        return (values=values, vectors=standardized_vectors)
    else
        if energies_window !== nothing && length(energies_window) == 2
            values = eigvals!(H, energies_window[1], energies_window[2])
        else
            values = eigvals!(H)
        end
        return (values=values, vectors=nothing)
    end
end

# Function to standardize eigenvector phases for consistency
function standardize_eigenvector_phases(vectors::AbstractMatrix{T}) where T
    standardized = copy(vectors)
    n_vecs = size(vectors, 2)
    
    for i in 1:n_vecs
        vec = @view standardized[:, i]
        
        # Find the index of the element with maximum absolute value
        max_idx = argmax(abs.(vec))
        max_element = vec[max_idx]
        
        if T <: Complex
            # For complex vectors, rotate phase to make max element real and positive
            if abs(max_element) > 1e-12
                phase_factor = conj(max_element) / abs(max_element)
                vec .*= phase_factor
            end
        else
            # For real vectors, ensure max element is positive
            if real(max_element) < 0
                vec .*= -1
            end
        end
    end
    
    return standardized
end

function convert_eigenvalues_to_harmonic_frequency_units(eigenvalues, harmonic_freq)
    """
    Convert eigenvalues to harmonic frequency units matching the time evolution.
    This ensures consistency between energy and time units in quench dynamics.
    """
    # Convert harmonic frequency to atomic units (Hartree)
    omega_au = Units.frequency_to_energy_au_value(harmonic_freq)
    
    # Convert eigenvalues (assumed to be in Hartree) to harmonic frequency units
    return eigenvalues ./ omega_au
end



end # module
