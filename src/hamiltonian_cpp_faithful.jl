module HamiltonianCpp

using Statistics
using LinearAlgebra
using SparseArrays
using Unitful # Add Unitful to the scope
using ..Units # Assuming Units module is in the same parent directory
using ..Utility # Add Utility module
using ..Parameters # To access PARAMS for Unitful quantities
using ..PhysicalConstantsJL # To access physical constants
using ..HamFunctionsCpp # Use functions from our C++ faithful implementation
using ..Profiling # Add profiling support

export buildHamiltonian, addRotation, addElectricField, addMagneticField, addSpinRotation, addElectricDipolarInteraction, addMagneticDipolarInteraction

function buildHamiltonian(channels::Vector{Vector{Float64}}, QN_max::Int, params::Parameters.Params, rotation::Bool=true; skip_terms::Set{Symbol}=Set{Symbol}())
    # Parse the hamiltonian dtype from params
    T = Utility.parse_hamiltonian_dtype(params.hamiltonian_dtype)

    set_matrix_elements_to_zero = false
    
    # Create the Hamiltonian matrix with specified type
    H = spzeros(T, QN_max, QN_max)
    
    # Track initial matrix memory and total system memory
    matrix_memory = Profiling.get_memory_size(H)
    total_memory = Base.gc_live_bytes()
    Profiling.track_memory("Hamiltonian_Matrix_Allocation", matrix_memory)
    Profiling.track_memory("Hamiltonian_Matrix_Allocation_Total", total_memory)
    
    # Get basic parameters using params
    basis_size = length(channels)
    
    println("Building Hamiltonian matrix ($(QN_max)×$(QN_max), type: $T)")
    println("Basis size: $basis_size channels")
    if !isempty(skip_terms)
        println("Skipping terms: $skip_terms")
    end
    
    # Add rotational terms
    if rotation && !(:rotation in skip_terms)
        Profiling.@time_operation "Hamiltonian_Rotation_Terms" begin
            addRotation(H, channels, params, set_matrix_elements_to_zero)
        end
        # Matrix size remains the same after in-place operations
        matrix_memory_after = Profiling.get_memory_size(H)
        total_memory_after = Base.gc_live_bytes()
        Profiling.track_memory("Hamiltonian_Rotation_Terms", matrix_memory_after)
        Profiling.track_memory("Hamiltonian_Rotation_Terms_Total", total_memory_after)
        println("Rotation matrix added to base hamiltonian.")
    end
    
    # Add electric field terms
    if ustrip(params.E) != 0 && !(:electric_field in skip_terms)
        Profiling.@time_operation "Hamiltonian_Electric_Field" begin
            addElectricField(H, channels, params, set_matrix_elements_to_zero, set_matrix_elements_to_zero)
        end
        matrix_memory_after = Profiling.get_memory_size(H)
        total_memory_after = Base.gc_live_bytes()
        Profiling.track_memory("Hamiltonian_Electric_Field", matrix_memory_after)
        Profiling.track_memory("Hamiltonian_Electric_Field_Total", total_memory_after)
        println("Electric matrix added to base hamiltonian.")
    end
    
    # Add magnetic field terms
    if ustrip(params.B) != 0 && !(:magnetic_field in skip_terms)
        Profiling.@time_operation "Hamiltonian_Magnetic_Field" begin
            addMagneticField(H, channels, params, set_matrix_elements_to_zero)
        end
        matrix_memory_after = Profiling.get_memory_size(H)
        total_memory_after = Base.gc_live_bytes()
        Profiling.track_memory("Hamiltonian_Magnetic_Field", matrix_memory_after)
        Profiling.track_memory("Hamiltonian_Magnetic_Field_Total", total_memory_after)
        println("Magnetic matrix added to base hamiltonian.")
    end
    
    # Add spin-rotation terms
    if ustrip(params.A) != 0 && !(:spin_rotation in skip_terms)
        Profiling.@time_operation "Hamiltonian_Spin_Rotation" begin
            addSpinRotation(H, channels, params, set_matrix_elements_to_zero)
        end
        matrix_memory_after = Profiling.get_memory_size(H)
        total_memory_after = Base.gc_live_bytes()
        Profiling.track_memory("Hamiltonian_Spin_Rotation", matrix_memory_after)
        Profiling.track_memory("Hamiltonian_Spin_Rotation_Total", total_memory_after)
        println("Spin-rotation matrix added to base hamiltonian.")
    end
    
    # Add electric dipolar interaction
    if ustrip(params.d_el) != 0 && length(params.interaction_geometries) > 0 && !(:electric_dipolar in skip_terms)
        Profiling.@time_operation "Hamiltonian_Electric_Dipolar" begin
            addElectricDipolarInteraction(H, channels, params, set_matrix_elements_to_zero)
        end
        matrix_memory_after = Profiling.get_memory_size(H)
        total_memory_after = Base.gc_live_bytes()
        Profiling.track_memory("Hamiltonian_Electric_Dipolar", matrix_memory_after)
        Profiling.track_memory("Hamiltonian_Electric_Dipolar_Total", total_memory_after)
        println("Electric dipolar matrix added to base hamiltonian.")
    end
    
    # Add magnetic dipolar interaction
    if ustrip(params.d_mg) != 0 && length(params.interaction_geometries) > 0 && !(:magnetic_dipolar in skip_terms)
        Profiling.@time_operation "Hamiltonian_Magnetic_Dipolar" begin
            addMagneticDipolarInteraction(H, channels, params, set_matrix_elements_to_zero)
        end
        matrix_memory_after = Profiling.get_memory_size(H)
        total_memory_after = Base.gc_live_bytes()
        Profiling.track_memory("Hamiltonian_Magnetic_Dipolar", matrix_memory_after)
        Profiling.track_memory("Hamiltonian_Magnetic_Dipolar_Total", total_memory_after)
        println("Magnetic dipolar matrix added to base hamiltonian.")
    end

    # Hermitize the Hamiltonian using the updated function from HamFunctionsCpp
    println("Starting Hermitizing the Hamiltonian matrix...")
    # Calculate the sparsity of the matrix

    nz_elements = nnz(H)
    sparsity = 1 - nz_elements / (QN_max * QN_max)
    println("Matrix sparsity: $sparsity (non-zero elements: $nz_elements out of $(QN_max * QN_max))")

    Profiling.@time_operation "Hamiltonian_Hermitization" begin
        HamFunctionsCpp.hermitize_matrix!(H)
    end
    matrix_memory_after = Profiling.get_memory_size(H)
    total_memory_after = Base.gc_live_bytes()
    Profiling.track_memory("Hamiltonian_Hermitization", matrix_memory_after)
    Profiling.track_memory("Hamiltonian_Hermitization_Total", total_memory_after)
    println("Hamiltonian matrix hermitized.")
    
    # Track final Hamiltonian memory
    final_matrix_memory = Profiling.get_memory_size(H)
    final_total_memory = Base.gc_live_bytes()
    Profiling.track_memory("Hamiltonian_Matrix_Final", final_matrix_memory)
    Profiling.track_memory("Hamiltonian_Matrix_Final_Total", final_total_memory)
    
    println("Hamiltonian construction complete!")

    return H
end

function addRotation(H::SparseMatrixCSC, channels::Vector{Vector{Float64}}, params::Parameters.Params, set_matrix_elements_to_zero::Bool=false, print_perm::Bool=true)
    if set_matrix_elements_to_zero
        H .= 0
    end

    if params.profile
        all_elements = []
    end
    
    basis_size = length(channels)
    if print_perm
        if params.print_colored
            printstyled("  Adding ")
            printstyled("rotational terms", color=:cyan)
            printstyled(" for ")
            printstyled("$basis_size", color=:cyan)
            printstyled(" basis states...\n")
        else    
            println("  Adding rotational terms for $basis_size states...")
        end
    end

    count = 0
    for i in 1:basis_size
        for molecule in 1:params.how_many_molecules
            j = channels[i][1 + params.QN_per_mol*(molecule-1)]
            H[i, i] += HamFunctionsCpp.rotationalFunction(j, params.Brot)

            if params.profile
                push!(all_elements, H[i, i])
            end
            count += 1
        end
    end

    if print_perm
        println("    Added $count rotational terms")
        if params.profile && !isempty(all_elements)
            if params.hamiltonian_dtype == "ComplexF64" || params.hamiltonian_dtype == "ComplexF32"
                all_elements = real.(all_elements)
            end
            mean_val = mean(all_elements)
            std_dev = std(all_elements)
            median_val = median(all_elements)
            min_val = minimum(all_elements)
            max_val = maximum(all_elements)
            println("    Rotational interaction terms stats")
            println("      Mean: $mean_val")
            println("      Std Dev: $std_dev")
            println("      Median: $median_val")
            println("      Min: $min_val")
            println("      Max: $max_val")
        end
    end
end

function addElectricField(H::SparseMatrixCSC, channels::Vector{Vector{Float64}}, params::Parameters.Params, set_matrix_elements_to_zero::Bool=false, print_perm::Bool=true)
    if ustrip(params.E) == 0 || ustrip(params.d_el) == 0
        return
    end
    
    if set_matrix_elements_to_zero
        H .= 0
    end

    if params.profile
        all_elements = []
    end
    
    basis_size = length(channels)
    if print_perm
        if params.print_colored
            printstyled("  Adding ")
            printstyled("electric field", color=:cyan)
            printstyled(" terms ")
            printstyled("(E = $(params.E))...\n", color=:cyan)
        else
            println("  Adding electric field terms (E = $(params.E))...")
        end
    end
    
    count = 0
    for i in 1:basis_size
        for iprim in 1:basis_size
            condition_for_j = HamFunctionsCpp.onlySinglejDifferByOne(channels, i, iprim, params.how_many_molecules, params.QN_per_mol)
            
            if condition_for_j[1]
                this_molecule_met_the_condition = condition_for_j[2]
                
                j = channels[i][1 + (this_molecule_met_the_condition-1) * params.QN_per_mol]
                jprim = channels[iprim][1 + (this_molecule_met_the_condition-1) * params.QN_per_mol]
                m = channels[iprim][2 + (this_molecule_met_the_condition-1) * params.QN_per_mol]
                
                H[i, iprim] += HamFunctionsCpp.electricFunction(params.d_el, params.E, j, jprim, m)

                if params.profile
                    push!(all_elements, H[i, iprim])
                end
                count += 1
            end
        end
    end
    
    if print_perm
        println("    Added $count electric field terms")
        if params.profile && !isempty(all_elements)
            if params.hamiltonian_dtype == "ComplexF64" || params.hamiltonian_dtype == "ComplexF32"
                all_elements = real.(all_elements)
            end
            mean_val = mean(all_elements)
            std_dev = std(all_elements)
            median_val = median(all_elements)
            min_val = minimum(all_elements)
            max_val = maximum(all_elements)
            println("    Electric field interaction terms stats")
            println("      Mean: $mean_val")
            println("      Std Dev: $std_dev")
            println("      Median: $median_val")
            println("      Min: $min_val")
            println("      Max: $max_val")
        end
    end
end

function addMagneticField(H::SparseMatrixCSC, channels::Vector{Vector{Float64}}, params::Parameters.Params, set_matrix_elements_to_zero::Bool=false, print_perm::Bool=true)
    if ustrip(params.B) == 0 || ustrip(params.d_mg) == 0
        return
    end
    
    if set_matrix_elements_to_zero
        H .= 0
    end

    if params.profile
        all_elements = []
    end
    
    basis_size = length(channels)
    if print_perm
        if params.print_colored
            printstyled("  Adding ")
            printstyled("magnetic field", color=:cyan)
            printstyled(" terms ")
            printstyled("(B = $(params.B))...\n", color=:cyan)
        else
            println("  Adding magnetic field terms (B = $(params.B))...")
        end
    end
    
    count = 0
    for i in 1:basis_size
        sum_of_mss = 0.0
        for molecule in 1:params.how_many_molecules
            ms = channels[i][4 + params.QN_per_mol*(molecule-1)]
            sum_of_mss += ms
        end
        H[i, i] += HamFunctionsCpp.magneticFunction(2u"μB", params.B, sum_of_mss)  # <-- Changed wrong implementation from d_mg to 2μB

        if params.profile
            push!(all_elements, H[i, i])
        end
        count += 1
    end
    
    if print_perm
        println("    Added $count magnetic field terms")
        if params.profile && !isempty(all_elements)
            if params.hamiltonian_dtype == "ComplexF64" || params.hamiltonian_dtype == "ComplexF32"
                all_elements = real.(all_elements)
            end
            mean_val = mean(all_elements)
            std_dev = std(all_elements)
            median_val = median(all_elements)
            min_val = minimum(all_elements)
            max_val = maximum(all_elements)
            println("    Magnetic field interaction terms stats")
            println("      Mean: $mean_val")
            println("      Std Dev: $std_dev")
            println("      Median: $median_val")
            println("      Min: $min_val")
            println("      Max: $max_val")
        end
    end
end

function addSpinRotation(H::SparseMatrixCSC, channels::Vector{Vector{Float64}}, params::Parameters.Params, set_matrix_elements_to_zero::Bool=false, print_perm::Bool=true)
    if ustrip(params.A) == 0
        return
    end
    
    if set_matrix_elements_to_zero
        H .= 0
    end
    
    basis_size = length(channels)
    spinrot_strength_au = Units.convert_param_to_au(params.A, :A_SR)

    if print_perm
        if params.print_colored
            printstyled("  Adding ")
            printstyled("spin-rotation", color=:cyan)
            printstyled(" terms ")
            printstyled("(A = $(params.A))...\n", color=:cyan)
        else
            println("  Adding spin-rotation terms (A = $(params.A))...")
        end
    end

    if params.profile
        all_elements = []
    end
    
    # Check if electric and magnetic field orientations are the same
    phi_el = ustrip(params.Phi_el)
    theta_el = ustrip(params.Theta_el)
    phi_mg = ustrip(params.Phi_mg)
    theta_mg = ustrip(params.Theta_mg)
    
    if phi_el == phi_mg && theta_el == theta_mg
        diagonal_count = 0
        offdiagonal_count = 0
        
        for i in 1:basis_size
            # Diagonal part
            H[i, i] += HamFunctionsCpp.spinRotationDiagonalFunction(params.A, channels, i, params.how_many_molecules)

            if params.profile
                push!(all_elements, H[i, i])
            end
            diagonal_count += 1
            
            for iprim in 1:basis_size
                if iprim != i
                    condition = HamFunctionsCpp.onlySinglemAndmsDiffer(channels, i, iprim, params.how_many_molecules, params.QN_per_mol)
                    
                    if condition[1]
                        this_molecule_met_the_condition = condition[2]
                        
                        m = channels[i][2 + (this_molecule_met_the_condition-1) * params.QN_per_mol]
                        mprim = channels[iprim][2 + (this_molecule_met_the_condition-1) * params.QN_per_mol]
                        ms = channels[i][4 + (this_molecule_met_the_condition-1) * params.QN_per_mol]
                        msprim = channels[iprim][4 + (this_molecule_met_the_condition-1) * params.QN_per_mol]
                        
                        if HamFunctionsCpp.msAndmssOffDiagonalSpinRotation(m, mprim, ms, msprim)
                            H[i, iprim] += spinrot_strength_au * 0.5

                            if params.profile
                                push!(all_elements, H[i, iprim])
                            end

                            offdiagonal_count += 1
                        end
                    end
                end
            end
        end
        
        if print_perm
            println("    Added $diagonal_count diagonal + $offdiagonal_count off-diagonal spin-rotation terms")
            if params.profile && !isempty(all_elements)
                if params.hamiltonian_dtype == "ComplexF64" || params.hamiltonian_dtype == "ComplexF32"
                    all_elements = real.(all_elements)
                end
                mean_val = mean(all_elements)
                std_dev = std(all_elements)
                median_val = median(all_elements)
                min_val = minimum(all_elements)
                max_val = maximum(all_elements)
                println("    Spin-rotation interaction terms stats")
                println("      Mean: $mean_val")
                println("      Std Dev: $std_dev")
                println("      Median: $median_val")
                println("      Min: $min_val")
                println("      Max: $max_val")
            end
        end
    else
        if print_perm
            println("    Different electric and magnetic field orientations not yet implemented in spin-rotation coupling")
        end
    end
end

# Coordinate transformation functions - faithful ports from C++ geometry.cpp
function polarAngleFromLaboratoryAngles(phi_R_quantity::Unitful.AbstractQuantity, theta_F_quantity::Unitful.AbstractQuantity, phi_F_quantity::Unitful.AbstractQuantity)
    # Convert all angles to radians using Unitful - this handles the unit conversion automatically
    phi_R_rad = Units.to_radians_value(phi_R_quantity)
    theta_F_rad = Units.to_radians_value(theta_F_quantity) 
    phi_F_rad = Units.to_radians_value(phi_F_quantity)
    
    # Apply the C++ formula: theta = acos(cos(phi_R) * sin(theta_F) * cos(phi_F) + sin(phi_R) * sin(theta_F) * sin(phi_F))
    theta_rad = acos(cos(phi_R_rad) * sin(theta_F_rad) * cos(phi_F_rad) + sin(phi_R_rad) * sin(theta_F_rad) * sin(phi_F_rad))
    
    return theta_rad  # Return in radians as Float64
end

function azimuthalAngleFromLaboratoryAngles(phi_R_quantity::Unitful.AbstractQuantity, theta_F_quantity::Unitful.AbstractQuantity, phi_F_quantity::Unitful.AbstractQuantity)
    # Convert all angles to radians using Unitful
    phi_R_rad = Units.to_radians_value(phi_R_quantity)
    theta_F_rad = Units.to_radians_value(theta_F_quantity)
    phi_F_rad = Units.to_radians_value(phi_F_quantity)
    
    # Apply the C++ formula: phi = -atan(cos(theta_F) * ((1 + tan(phi_R) * tan(phi_F)) / (tan(phi_F) - tan(phi_R))))
    numerator = 1 + tan(phi_R_rad) * tan(phi_F_rad)
    denominator = tan(phi_F_rad) - tan(phi_R_rad)
    
    # Handle potential division by zero
    if abs(denominator) < 1e-12
        return 0.0  # Default fallback
    end
    
    phi_rad = -atan(cos(theta_F_rad) * (numerator / denominator))
    
    return phi_rad  # Return in radians as Float64
end

function addElectricDipolarInteraction(H::SparseMatrixCSC, channels::Vector{Vector{Float64}}, params::Parameters.Params, set_matrix_elements_to_zero::Bool=false, print_perm::Bool=true)
    if ustrip(params.d_el) == 0
        return
    end

    if params.profile
        all_elements = []
    end
    
    if set_matrix_elements_to_zero
        H .= 0
    end
    
    basis_size = length(channels)

    if print_perm
        if params.print_colored
            printstyled("  Adding ")
            printstyled("electric dipolar", color=:cyan)
            printstyled(" interactions ")
            printstyled("(d_el = $(params.d_el))...\n", color=:cyan)
        else
            println("  Adding electric dipolar interactions (d_el = $(params.d_el))...")
            println("    Interaction geometries:")
            println("       R_1 = 0 nm, φ_R1 = 0°")
            for (num, (R, phi_R)) in enumerate(params.interaction_geometries)
                println("       R_$(num+1) = $(ustrip(R)) nm, φ_R$(num+1) = $(ustrip(phi_R))°")
            end
        end
    end
    
    count = 0
    for i in 1:basis_size
        for iprim in 1:basis_size
            condition_for_j_m_pair = HamFunctionsCpp.onlySinglePairOfjmMayDiffer(channels, i, iprim, params.how_many_molecules, params.QN_per_mol)
            
            if condition_for_j_m_pair[1]
                this_pair_met_the_condition = condition_for_j_m_pair[2]
                molecule1 = this_pair_met_the_condition[1]
                molecule2 = this_pair_met_the_condition[2]
                
                # Extract quantum numbers for molecule1
                j1_idx = (molecule1-1) * params.QN_per_mol + 1
                m1_idx = (molecule1-1) * params.QN_per_mol + 2
                j1 = channels[i][j1_idx]
                m1 = channels[i][m1_idx]
                j1prim = channels[iprim][j1_idx]
                m1prim = channels[iprim][m1_idx]
                
                # Extract quantum numbers for molecule2
                j2_idx = (molecule2-1) * params.QN_per_mol + 1
                m2_idx = (molecule2-1) * params.QN_per_mol + 2
                j2 = channels[i][j2_idx]
                m2 = channels[i][m2_idx]
                j2prim = channels[iprim][j2_idx]
                m2prim = channels[iprim][m2_idx]
                
                # Validate quantum numbers
                if abs(m1) > j1 || abs(m1prim) > j1prim || abs(m2) > j2 || abs(m2prim) > j2prim
                    continue
                end
                
                # Get relative distance and angle from interaction_geometries
                if molecule2 <= length(params.interaction_geometries) + 1
                    pair_idx = molecule2 - 1
                    if pair_idx >= 1 && pair_idx <= length(params.interaction_geometries)
                        R_quantity, relative_lab_angle_quantity = params.interaction_geometries[pair_idx]
                        
                        theta = polarAngleFromLaboratoryAngles(relative_lab_angle_quantity, params.Theta_el, params.Phi_el)
                        phi = azimuthalAngleFromLaboratoryAngles(relative_lab_angle_quantity, params.Theta_el, params.Phi_el)
                        
                        interaction_term = HamFunctionsCpp.electricDipolarFunction(
                            params.d_el, R_quantity, theta, phi,
                            j1, j1prim, m1, m1prim,
                            j2, j2prim, m2, m2prim
                        )
                        
                        H[i, iprim] += interaction_term
                        
                        if params.profile
                            push!(all_elements, H[i, iprim])
                        end
                        count += 1
                    end
                end
            end
        end
    end
    
    if print_perm
        println("    Added $count electric dipolar interaction terms")
        if params.profile && !isempty(all_elements)
            if params.hamiltonian_dtype == "ComplexF64" || params.hamiltonian_dtype == "ComplexF32"
                all_elements = real.(all_elements)
            end
            mean_val = mean(all_elements)
            std_dev = std(all_elements)
            median_val = median(all_elements)
            min_val = minimum(all_elements)
            max_val = maximum(all_elements)
            println("    Electric dipolar interaction terms stats")
            println("      Mean: $mean_val")
            println("      Std Dev: $std_dev")
            println("      Median: $median_val")
            println("      Min: $min_val")
            println("      Max: $max_val")
        end
    end
end

function addMagneticDipolarInteraction(H::SparseMatrixCSC, channels::Vector{Vector{Float64}}, params::Parameters.Params, set_matrix_elements_to_zero::Bool=false, print_perm::Bool=true)
    if ustrip(params.d_mg) == 0
        return
    end
    
    if set_matrix_elements_to_zero
        H .= 0
    end
    
    basis_size = length(channels)

    if print_perm
        if params.print_colored
            printstyled("  Adding ")
            printstyled("magnetic dipolar", color=:cyan)
            printstyled(" interactions ")
            printstyled("(d_mg = $(ustrip(params.d_mg / 1u"μB")) μB)...\n", color=:cyan)
        else
            println("  Adding magnetic dipolar interactions (d_mg = $(ustrip(params.d_mg / 1u"μB")) μB)...")
            println("    Interaction geometries:")
            println("       R_1 = 0 nm, φ_R1 = 0°")
            for (num, (R, phi_R)) in enumerate(params.interaction_geometries)
                println("       R_$(num+1) = $(ustrip(R)) nm, φ_R$(num+1) = $(ustrip(phi_R))°")
            end
        end
    end

    if params.profile
        all_elements = []
    end
    
    count = 0
    for i in 1:basis_size
        for iprim in 1:basis_size
            condition_for_ms_pair = HamFunctionsCpp.onlySinglePairOfmsMayDifferMaxByOne(channels, i, iprim, params.how_many_molecules, params.QN_per_mol)
            
            if condition_for_ms_pair[1]
                this_pair_met_the_condition = condition_for_ms_pair[2]
                molecule1 = this_pair_met_the_condition[1]
                molecule2 = this_pair_met_the_condition[2]
                
                # Extract ms quantum numbers
                ms1 = channels[i][4 + (molecule1-1) * params.QN_per_mol]
                ms1prim = channels[iprim][4 + (molecule1-1) * params.QN_per_mol]
                ms2 = channels[i][4 + (molecule2-1) * params.QN_per_mol]
                ms2prim = channels[iprim][4 + (molecule2-1) * params.QN_per_mol]
                
                # Get spin values
                s1 = channels[i][3 + (molecule1-1) * params.QN_per_mol]
                s2 = channels[i][3 + (molecule2-1) * params.QN_per_mol]
                
                # Get relative distance and angle from interaction_geometries
                if molecule2 <= length(params.interaction_geometries) + 1
                    pair_idx = molecule2 - 1
                    if pair_idx >= 1 && pair_idx <= length(params.interaction_geometries)
                        R_quantity, relative_lab_angle_quantity = params.interaction_geometries[pair_idx]
                        
                        theta = polarAngleFromLaboratoryAngles(relative_lab_angle_quantity, params.Theta_mg, params.Phi_mg)
                        phi = azimuthalAngleFromLaboratoryAngles(relative_lab_angle_quantity, params.Theta_mg, params.Phi_mg)
                        
                        interaction_term = HamFunctionsCpp.magneticDipolarFunction(
                            params.d_mg, R_quantity, theta, phi,
                            s1, ms1, ms1prim,
                            s2, ms2, ms2prim
                        )
                        
                        H[i, iprim] += interaction_term
                        
                        if params.profile
                            push!(all_elements, H[i, iprim])
                        end
                        count += 1
                    end
                end
            end
        end
    end
    
    if print_perm
        println("    Added $count magnetic dipolar interaction terms")
        if params.profile && !isempty(all_elements)
            if params.hamiltonian_dtype == "ComplexF64" || params.hamiltonian_dtype == "ComplexF32"
                all_elements = real.(all_elements)
            end
            mean_val = mean(all_elements)
            std_dev = std(all_elements)
            median_val = median(all_elements)
            min_val = minimum(all_elements)
            max_val = maximum(all_elements)
            println("    Magnetic dipolar interaction terms stats")
            println("      Mean: $mean_val")
            println("      Std Dev: $std_dev")
            println("      Median: $median_val")
            println("      Min: $min_val")
            println("      Max: $max_val")
        end
    end
end

# Function to create a matrix containing only a specific term
function createTermMatrix(channels::Vector{Vector{Float64}}, QN_max::Int, params::Parameters.Params, term::Symbol, print_perm=false)
    T = Utility.parse_hamiltonian_dtype(params.hamiltonian_dtype)
    matrix = spzeros(T, QN_max, QN_max)
    
    if term == :electric_field
        addElectricField(matrix, channels, params, true, print_perm)
    elseif term == :magnetic_field
        addMagneticField(matrix, channels, params, true, print_perm)
    elseif term == :spin_rotation
        addSpinRotation(matrix, channels, params, true, print_perm)
    elseif term == :electric_dipolar
        addElectricDipolarInteraction(matrix, channels, params, true, print_perm)
    elseif term == :magnetic_dipolar
        addMagneticDipolarInteraction(matrix, channels, params, true, print_perm)
    else
        error("Unknown term: $term")
    end
    
    return matrix
end

# Function to create a matrix containing multiple terms
function createCombinedMatrix(channels::Vector{Vector{Float64}}, QN_max::Int, params::Parameters.Params, terms::Vector{Symbol})
    T = Utility.parse_hamiltonian_dtype(params.hamiltonian_dtype)
    matrix = spzeros(T, QN_max, QN_max)
    
    for term in terms
        if term == :electric_field
            addElectricField(matrix, channels, params, false, false)
        elseif term == :magnetic_field
            addMagneticField(matrix, channels, params, false, false)
        elseif term == :spin_rotation
            addSpinRotation(matrix, channels, params, false, false)
        elseif term == :electric_dipolar
            addElectricDipolarInteraction(matrix, channels, params, false, false)
        elseif term == :magnetic_dipolar
            addMagneticDipolarInteraction(matrix, channels, params, false, false)
        else
            @warn "Unknown term: $term"
        end
    end
    
    return matrix
end

end # module HamiltonianCpp