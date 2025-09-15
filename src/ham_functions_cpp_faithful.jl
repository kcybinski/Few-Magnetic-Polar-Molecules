module HamFunctionsCpp

using LinearAlgebra
using SparseArrays
using Unitful # Add Unitful to the scope
using ..Units # Assuming Units module is in the same parent directory
using ..Parameters # To access PARAMS for Unitful quantities
using ..PhysicalConstantsJL # To access physical constants
using WignerSymbols # Make WignerSymbols.jl a hard dependency

# Helper function to convert C++ function calls to Julia equivalents
function ClebschGordanCoefficient(j1::Int, m1::Int, j2::Int, m2::Int, Jtot::Int, Mtot::Int)
    return WignerSymbols.clebschgordan(j1, m1, j2, m2, Jtot, Mtot)
end

function electricDipoleOperator(q::Int, j::Float64, jprim::Float64, m::Float64, mprim::Float64)
    # without d_el
    result = 0.0
    
    # Validate quantum numbers before proceeding
    if abs(m) > j || abs(mprim) > jprim
        return 0.0  # Invalid quantum numbers: |m| > j or |mprim| > jprim
    end
    
    # Ensure m_prim = m + q for CG(j,m,1,q, j_prim, m_prim) to be non-zero
    if mprim != m + q
        return 0.0
    end
    
    # Ensure j_prim is in [j-1, j, j+1] for CGs to be non-zero
    if abs(j - jprim) > 1
        return 0.0
    end
    
    # Convert to Int for Clebsch-Gordan coefficients only when needed
    j_int = Int(j)
    jprim_int = Int(jprim)
    m_int = Int(m)
    mprim_int = Int(mprim)
    
    result = ClebschGordanCoefficient(j_int, m_int, 1, q, jprim_int, mprim_int) * 
             ClebschGordanCoefficient(j_int, 0, 1, 0, jprim_int, 0) * 
             sqrt((2 * j + 1) / (2 * jprim + 1))  # Use Float64 directly for arithmetic
    
    return result
end

function spinLadderOperator(q::Int, s::Float64, ms::Float64, msprim::Float64)
    # without d_mg
    result = 0.0
    
    if q == 0
        result = ms
    elseif abs(q) == 1
        result = -1.0 / sqrt(2.0) * sqrt(s * (s + 1) - ms * msprim)
    else
        println("ERROR! q in spinLadder Operator is larger than 1.")
    end
    
    return result
end

function rotationalFunction(j::Float64, B_rot_quantity::Unitful.AbstractQuantity)
    B_rot_au = Units.convert_param_to_au(B_rot_quantity, :Brot)
    result = B_rot_au * j * (j + 1)
    return result
end

function electricFunction(d_el_quantity::Unitful.AbstractQuantity, el_field_quantity::Unitful.AbstractQuantity, j::Float64, jprim::Float64, m::Float64)
    d_el_au = Units.convert_param_to_au(d_el_quantity, :d_el)
    el_field_au = Units.convert_param_to_au(el_field_quantity, :E_field)
    result = -d_el_au * el_field_au * electricDipoleOperator(0, j, jprim, m, m)
    return result
end

function magneticFunction(d_mg_quantity::Unitful.AbstractQuantity, mg_field_quantity::Unitful.AbstractQuantity, sum_of_mss::Float64)
    d_mg_au = Units.convert_param_to_au(d_mg_quantity, :d_mg)
    mg_field_au = Units.convert_param_to_au(mg_field_quantity, :B_field)
    return d_mg_au * sum_of_mss * mg_field_au
end

function spinRotationDiagonalFunction(spinrot_strength_quantity::Unitful.AbstractQuantity, channels::Vector{Vector{Float64}}, i::Int, how_many_molecules::Int)
    spinrot_strength_au = Units.convert_param_to_au(spinrot_strength_quantity, :A_SR)
    sum = 0.0
    
    for molecule in 1:how_many_molecules
        m = channels[i][2 + 4*(molecule-1)]
        ms = channels[i][4 + 4*(molecule-1)]
        sum += m * ms
    end
    
    return spinrot_strength_au * sum
end

function electricDipolarFunction(d_el_quantity::Unitful.AbstractQuantity, R_quantity::Unitful.AbstractQuantity, theta::Float64, phi::Float64, 
                                j1::Float64, j1prim::Float64, m1::Float64, m1prim::Float64, 
                                j2::Float64, j2prim::Float64, m2::Float64, m2prim::Float64)
    d_el_au = Units.convert_param_to_au(d_el_quantity, :d_el)
    R_au = Units.convert_param_to_au(R_quantity, :R_dist)
    
    dq1dq2 = 0.0
    q1 = Int(m1prim - m1)  # Only convert to Int when needed for q calculation
    q2 = Int(m2prim - m2)
    angular_part = 0.0 + 0.0im
    
    # Add validation: electric dipole transitions require |q| ≤ 1
    if abs(q1) > 1 || abs(q2) > 1
        return 0.0 + 0.0im  # Invalid electric dipole transition
    end
    
    dq1dq2 = electricDipoleOperator(q1, j1, j1prim, m1, m1prim) * electricDipoleOperator(q2, j2, j2prim, m2, m2prim)
    
    if dq1dq2 == 0.0
        return 0.0
    end
    
    angular_part = 0.0 + 0.0im
    
    if q1 == 0 && q2 == 0
        angular_part = 1 - 3 * cos(theta)^2
    elseif (q1 == -1 && q2 == 1) || (q1 == 1 && q2 == -1)
        angular_part = (1 - 3 * cos(theta)^2) / 2
    elseif (q1 == 0 && q2 == 1) || (q1 == 1 && q2 == 0)
        angular_part = 3 / sqrt(2) * sin(theta) * cos(theta) * exp(-im * phi)
    elseif (q1 == 0 && q2 == -1) || (q1 == -1 && q2 == 0)
        angular_part = 3 / sqrt(2) * sin(theta) * cos(theta) * exp(im * phi)
    elseif q1 == 1 && q2 == 1
        angular_part = -3 / 2 * sin(theta)^2 * exp(-im * (2 * phi))
    elseif q1 == -1 && q2 == -1
        angular_part = -3 / 2 * sin(theta)^2 * exp(im * (2 * phi))
    else
        return 0.0
    end
    
    return (d_el_au * d_el_au) / (R_au^3) * angular_part * dq1dq2
end

function magneticDipolarFunction(d_mg_quantity::Unitful.AbstractQuantity, R_quantity::Unitful.AbstractQuantity, theta::Float64, phi::Float64, 
                               s1::Float64, ms1::Float64, ms1prim::Float64, 
                               s2::Float64, ms2::Float64, ms2prim::Float64)
    d_mg_au = Units.convert_param_to_au(d_mg_quantity, :d_mg)
    R_au = Units.convert_param_to_au(R_quantity, :R_dist)
    
    sq1sq2 = 0.0
    q1 = Int(ms1prim - ms1)
    q2 = Int(ms2prim - ms2)
    angular_part = 0.0 + 0.0im
    
    sq1sq2 = spinLadderOperator(q1, s1, ms1, ms1prim) * spinLadderOperator(q2, s2, ms2, ms2prim)
    
    angular_part = 0.0 + 0.0im
    
    if q1 == 0 && q2 == 0
        angular_part = 1 - 3 * cos(theta)^2
    elseif (q1 == -1 && q2 == 1) || (q1 == 1 && q2 == -1)
        angular_part = (1 - 3 * cos(theta)^2) / 2
    elseif (q1 == 0 && q2 == 1) || (q1 == 1 && q2 == 0)
        angular_part = 3 / sqrt(2) * sin(theta) * cos(theta) * exp(-im * phi)
    elseif (q1 == 0 && q2 == -1) || (q1 == -1 && q2 == 0)
        angular_part = 3 / sqrt(2) * sin(theta) * cos(theta) * exp(im * phi)
    elseif q1 == 1 && q2 == 1
        angular_part = -3 / 2 * sin(theta)^2 * exp(-im * (2 * phi))
    elseif q1 == -1 && q2 == -1
        angular_part = -3 / 2 * sin(theta)^2 * exp(im * (2 * phi))
    else
        return 0.0
    end
    
    alpha_squared = PhysicalConstantsJL.fine_structure_constant^2
    
    return (alpha_squared * d_mg_au^2 / R_au^3) * angular_part * sq1sq2
end

function onlySinglejDifferByOne(channels::Vector{Vector{Float64}}, i::Int, iprim::Int, how_many_molecules::Int, QN_per_mol::Int)
    # j1 m1 S1 ms1 | j2 m2 S2 ms2 | j3 m3 S3 ms3 | j4 m4 S4 ms4
    
    everythingTheSameBesideSinglej = false
    this_molecule_met_the_condition = -1
    
    for molecule in 1:how_many_molecules
        # Check if everything before j for this molecule is the same
        uptoj = true
        start_idx = 1 + (molecule-1) * QN_per_mol
        for idx in 1:(start_idx-1)
            if channels[i][idx] != channels[iprim][idx]
                uptoj = false
                break
            end
        end
        
        # Check if everything after j for this molecule is the same
        afterj = true
        after_j_start = start_idx + 1  # Skip j, start from m
        for idx in after_j_start:length(channels[i])
            if channels[i][idx] != channels[iprim][idx]
                afterj = false
                break
            end
        end
        
        if uptoj && afterj
            everythingTheSameBesideSinglej = true
            this_molecule_met_the_condition = molecule
            break  # this condition can happen only once, for one molecule
        end
    end
    
    if !everythingTheSameBesideSinglej
        return (false, 10)
    else
        j_idx = 1 + (this_molecule_met_the_condition-1) * QN_per_mol
        j = channels[i][j_idx]      # Keep as Float64
        jprim = channels[iprim][j_idx]  # Keep as Float64
        
        if abs(j - jprim) == 1.0   # Compare as Float64
            return (true, this_molecule_met_the_condition)
        else
            return (false, 10)
        end
    end
end

function onlySinglemAndmsDiffer(channels::Vector{Vector{Float64}}, i::Int, iprim::Int, how_many_molecules::Int, QN_per_mol::Int)
    # j1 m1 S1 ms1 | j2 m2 S2 ms2 | j3 m3 S3 ms3 | j4 m4 S4 ms4
    
    everythingTheSameBesideSinglemAndms = false
    this_molecule_met_the_condition = 10
    
    for molecule in 1:how_many_molecules
        # Check if everything before m for this molecule is the same
        uptom = true
        start_idx = 1 + (molecule-1) * QN_per_mol
        m_idx = start_idx + 1  # m is at position 2 relative to j
        for idx in 1:(m_idx-1)
            if channels[i][idx] != channels[iprim][idx]
                uptom = false
                break
            end
        end
        
        # Check if everything after ms for this molecule is the same  
        afterms = true
        ms_idx = start_idx + 3  # ms is at position 4 relative to j
        for idx in (ms_idx+1):length(channels[i])
            if channels[i][idx] != channels[iprim][idx]
                afterms = false
                break
            end
        end
        
        if uptom && afterms
            everythingTheSameBesideSinglemAndms = true
            this_molecule_met_the_condition = molecule
            break  # this condition can happen only for one molecule
        end
    end
    
    return (everythingTheSameBesideSinglemAndms, this_molecule_met_the_condition)
end

function onlySinglePairOfjmMayDiffer(channels::Vector{Vector{Float64}}, i::Int, iprim::Int, how_many_molecules::Int, QN_per_mol::Int)
    pairs = Tuple{Int, Int}[]
    
    if how_many_molecules == 2
        push!(pairs, (1, 2))  # Julia 1-based indexing: (0,1) -> (1,2)
    elseif how_many_molecules == 3
        push!(pairs, (1, 2))
        push!(pairs, (1, 3))
        push!(pairs, (2, 3))
    elseif how_many_molecules == 4
        push!(pairs, (1, 2))
        push!(pairs, (1, 3))
        push!(pairs, (1, 4))
        push!(pairs, (2, 3))
        push!(pairs, (2, 4))
        push!(pairs, (3, 4))
    else
        println("Error! Unknown pairs for electric dipolar interaction.")
    end
    
    # Julia channel format: [j, m, S_actual, ms_actual] for each molecule
    for (mol1, mol2) in pairs
        # Calculate indices for j quantum numbers (position 1 for each molecule)
        j1_idx = (mol1-1) * QN_per_mol + 1  # j position for molecule 1
        j2_idx = (mol2-1) * QN_per_mol + 1  # j position for molecule 2
        
        j1 = channels[i][j1_idx]      # Keep as Float64
        j1prim = channels[iprim][j1_idx]  # Keep as Float64
        j2 = channels[i][j2_idx]      # Keep as Float64
        j2prim = channels[iprim][j2_idx]  # Keep as Float64
        
        # Get m quantum numbers for validation (position 2 for each molecule)
        m1_idx = (mol1-1) * QN_per_mol + 2  # m position for molecule 1
        m2_idx = (mol2-1) * QN_per_mol + 2  # m position for molecule 2
        m1 = channels[i][m1_idx]      # m1 
        m1prim = channels[iprim][m1_idx]  # m1prim 
        m2 = channels[i][m2_idx]      # m2 
        m2prim = channels[iprim][m2_idx]  # m2prim 
        
        # Validate quantum numbers
        if abs(m1) > j1 || abs(m1prim) > j1prim || abs(m2) > j2 || abs(m2prim) > j2prim
            return (false, (10, 10))  # Return false for invalid combinations
        end
        
        # Electric dipole selection rules: |Δj| = 1 AND |Δm| ≤ 1 for each molecule
        if abs(j1 - j1prim) == 1.0 && abs(j2 - j2prim) == 1.0  # Compare as Float64
            # Additional check: electric dipole selection rule |Δm| ≤ 1
            if abs(m1 - m1prim) <= 1.0 && abs(m2 - m2prim) <= 1.0
                # Check if everything before j1,m1 is the same (everything before molecule 1)
                uptoj1m1 = true
                for idx in 1:(j1_idx-1)
                    if channels[i][idx] != channels[iprim][idx]
                        uptoj1m1 = false
                        break
                    end
                end
                
                # Check if everything between j1,m1 and j2,m2 is the same (S1, ms1 and everything up to molecule 2)
                betweenj1m1andj2m2 = true
                start_between = m1_idx + 1  # After m1 (start from S1)
                end_between = j2_idx - 1   # Up to but not including j2
                for idx in start_between:end_between
                    if channels[i][idx] != channels[iprim][idx]
                        betweenj1m1andj2m2 = false
                        break
                    end
                end
                
                # Check if everything after j2,m2 is the same (S2, ms2 and any remaining molecules)
                afterj2m2 = true
                start_after = m2_idx + 1  # After m2 (start from S2)
                for idx in start_after:length(channels[i])
                    if channels[i][idx] != channels[iprim][idx]
                        afterj2m2 = false
                        break
                    end
                end
                
                if uptoj1m1 && betweenj1m1andj2m2 && afterj2m2
                    return (true, (mol1, mol2))  # these conditions can happen only once
                end
            end
        end
    end
    
    return (false, (10, 10))
end

function onlySinglePairOfmsMayDifferMaxByOne(channels::Vector{Vector{Float64}}, i::Int, iprim::Int, how_many_molecules::Int, QN_per_mol::Int)
    pairs = Tuple{Int, Int}[]
    
    if how_many_molecules == 2
        push!(pairs, (1, 2))  # Julia 1-based indexing: (0,1) -> (1,2)
    elseif how_many_molecules == 3
        push!(pairs, (1, 2))
        push!(pairs, (1, 3))
        push!(pairs, (2, 3))
    elseif how_many_molecules == 4
        push!(pairs, (1, 2))
        push!(pairs, (1, 3))
        push!(pairs, (1, 4))
        push!(pairs, (2, 3))
        push!(pairs, (2, 4))
        push!(pairs, (3, 4))
    else
        println("Error! Unknown pairs for magnetic dipolar interaction.")
    end
    
    # j1 m1 S1 ms1 | j2 m2 S2 ms2 | j3 m3 S3 ms3 | j4 m4 S4 ms4
    for (mol1, mol2) in pairs
        index_ms1 = QN_per_mol * (mol1-1) + 4  # ms position for molecule 1
        index_ms2 = QN_per_mol * (mol2-1) + 4  # ms position for molecule 2
        
        ms1 = channels[i][index_ms1]  # Keep as Float64 since ms can be half-integer
        ms1prim = channels[iprim][index_ms1]
        ms2 = channels[i][index_ms2]
        ms2prim = channels[iprim][index_ms2]
        
        # In C++, this checks abs(ms1 - ms1prim) == 2 because ms values are doubled
        # In Julia, ms values are actual, so we check abs(ms1 - ms1prim) == 1.0
        if (abs(ms1 - ms1prim) == 1.0 && abs(ms2 - ms2prim) == 1.0) || (ms1 == ms1prim && ms2 == ms2prim)
            # Check if everything before ms1 is the same
            uptoms1 = true
            for idx in 1:(index_ms1-1)
                if channels[i][idx] != channels[iprim][idx]
                    uptoms1 = false
                    break
                end
            end
            
            # Check if everything between ms1 and ms2 is the same (ms can differ)
            betweenms1andms2 = true
            start_between = index_ms1 + 1  # After ms1
            end_between = index_ms2 - 1   # Up to but not including ms2
            for idx in start_between:end_between
                if channels[i][idx] != channels[iprim][idx]
                    betweenms1andms2 = false
                    break
                end
            end
            
            # Check if everything after ms2 is the same (ms can differ)
            afterms2 = true
            start_after = index_ms2 + 1  # After ms2
            for idx in start_after:length(channels[i])
                if channels[i][idx] != channels[iprim][idx]
                    afterms2 = false
                    break
                end
            end
            
            if uptoms1 && betweenms1andms2 && afterms2
                return (true, (mol1, mol2))  # these conditions can happen only once
            end
        end
    end
    
    return (false, (10, 10))
end

function msAndmssOffDiagonalSpinRotation(m::Float64, mprim::Float64, ms::Float64, msprim::Float64)
    # for spinRotation off-diagonal terms
    if ms == msprim + 1.0 && m == mprim - 1.0
        return true
    elseif ms == msprim - 1.0 && m == mprim + 1.0
        return true
    else 
        return false
    end
end

function hermitize_matrix!(H::Matrix{T}) where T
    """
    Efficiently enforce Hermiticity by copying upper triangle to lower triangle.
    Modifies the matrix in-place for memory efficiency.
    """
    n = size(H, 1)
    
    if size(H, 2) != n
        error("Matrix must be square for Hermitization")
    end
    
    # Make diagonal elements real if complex type
    for i in 1:n
        if T <: Complex
            H[i, i] = real(H[i, i])
        end
    end
    
    # Copy upper triangle to lower triangle
    for i in 1:n
        for j in (i+1):n
            if T <: Complex
                H[j, i] = conj(H[i, j])
            else
                H[j, i] = H[i, j]  # Real symmetric matrix
            end
        end
    end
    
    return H
end

# Efficient hermitize for sparse matrix using only non-zero elements
function hermitize_matrix!(H::SparseMatrixCSC)
    n = size(H, 1)
    if size(H, 2) != n
        error("Matrix must be square for Hermitization")
    end

    T = eltype(H)
    
    # Create a copy of the upper triangular part (including diagonal)
    # This prevents modifying the matrix while iterating over it
    upper_indices = Tuple{Int,Int}[]
    upper_values = Vector{T}()
    
    # First pass: collect upper triangle elements and fix diagonal
    for col in 1:n
        for j in SparseArrays.nzrange(H, col)
            row = H.rowval[j]
            val = H.nzval[j]
            
            if row == col  # Diagonal element
                if T <: Complex
                    # Make diagonal elements real
                    H.nzval[j] = real(val)
                    push!(upper_indices, (row, col))
                    push!(upper_values, real(val))
                else
                    push!(upper_indices, (row, col))
                    push!(upper_values, val)
                end
            elseif row < col  # Upper triangle
                push!(upper_indices, (row, col))
                push!(upper_values, val)
            end
        end
    end
    
    # Second pass: fill lower triangle based on collected upper elements
    for (idx, (row, col)) in enumerate(upper_indices)
        if row != col  # Skip diagonal (already fixed)
            val = upper_values[idx]
            if T <: Complex
                H[col, row] = conj(val)
            else
                H[col, row] = val
            end
        end
    end
    
    return H
end

end