module BasisCppFaithful

"""
    generate_channels_cpp_faithful(Mtot_target::Int, print::Bool=false; JMAX::Int=1, S_input::Int=1, how_many_molecules::Int=1)

Faithful reproduction of the C++ generateChannels algorithm.

# Arguments
- `Mtot_target::Int`: Target total magnetic quantum number (as integer, matching C++ convention)
- `print::Bool=false`: Whether to print debug output
- `JMAX::Int=1`: Maximum rotational quantum number
- `S_input::Int=1`: Spin quantum number in C++ convention (2*actual_spin as integer)
- `how_many_molecules::Int=1`: Number of molecules in the system

# Returns
- `Vector{Int}`: Flat vector of quantum numbers [j1,m1,S1,ms1,j2,m2,S2,ms2,...] as integers

# Notes
The C++ convention stores:
- Spin S as 2*S (e.g., spin-1/2 → S=1, spin-1 → S=2)
- Magnetic spin ms as 2*ms (e.g., ms=±0.5 → ms=±1, ms=±1 → ms=±2)
- Mtot as integer (e.g., Mtot=0.5 → Mtot_target=1 when doubled)
"""

using ..Utility

export generate_channels_cpp_faithful

function check_mtot_constraint(m_tab, ms_tab, target_mtot_doubled; mode="simple sum", override=false)
    if override
        return true
    end

    if mode == "cpp"
        # This is here for legacy reasons, to reproduce the wrong C++ behavior
        # This seems to only work if Mtot is an integer
        ms_sum = sum(ms_tab)
        ms_sum_half = ms_sum * 0.5
        cpp_round_ms_sum = ms_sum_half >= 0 ? Int(floor(ms_sum_half + 0.5)) : Int(ceil(ms_sum_half - 0.5))
        mtot_calculated = sum(m_tab) + cpp_round_ms_sum
        mtot_calculated_doubled = round(Int, 2*mtot_calculated)
        return mtot_calculated_doubled == target_mtot_doubled
    elseif mode == "simple sum"
        ms_sum = sum(ms_tab) * 0.5
        m_sum = sum(m_tab)
        mtot_calculated = m_sum + ms_sum
        mtot_calculated_doubled = round(Int, 2 * mtot_calculated)
        return mtot_calculated_doubled == target_mtot_doubled
    else
        error("Unknown mode for check_mtot_constraint")
    end
end

function validate_Mtot_possibility(Mtot_target_doubled, JMax, Smax, how_many_molecules)
    # Calculate min and max possible Mtot values, considering all molecules
    Mtot_target = Mtot_target_doubled * 0.5  # Convert back to actual Mtot
    max_ms_per_mol = Smax * 0.5
    min_ms_per_mol = -Smax * 0.5
    min_mj_per_mol = -JMax
    max_mj_per_mol = JMax

    # Check if Mtot is within the accessible range
    min_mtot = min_mj_per_mol * how_many_molecules + min_ms_per_mol * how_many_molecules
    max_mtot = max_mj_per_mol * how_many_molecules + max_ms_per_mol * how_many_molecules

    if Mtot_target < min_mtot || Mtot_target > max_mtot
        error("Mtot $(Mtot_target) is outside accesible range! The minimal Mtot=$(min_mtot), maximal=$(max_mtot)")
    end

    integer_mtot_allowed = mod(max_ms_per_mol * how_many_molecules, 1) == 0 ? true : false
    half_integer_mtot_allowed = !integer_mtot_allowed

    allowed_mtot_values = integer_mtot_allowed ? collect(min_mtot:1:max_mtot) : collect(min_mtot:1:max_mtot)


    if integer_mtot_allowed && mod(Mtot_target, 1) != 0
        error("Mtot $(Mtot_target) is half-integer but only integer Mtot is allowed for spin=$(Smax)/2, Jmax=$(JMax) and $(how_many_molecules) molecules!\nAllowed Mtot values: $(allowed_mtot_values)")
    elseif half_integer_mtot_allowed && mod(Mtot_target, 1) == 0
        error("Mtot $(Mtot_target) is integer but only half-integer Mtot is allowed for spin=$(Smax)/2, Jmax=$(JMax) and $(how_many_molecules) molecules!\nAllowed Mtot values: $(allowed_mtot_values)")
    end
    return true
end

    

function generate_channels_cpp_faithful(Mtot_target::Int, print_debug::Bool=false; JMAX::Int=1, S_input::Int=1, how_many_molecules::Int=1, QN_per_mol::Int=4, constraint_check_mode="simple sum", override_mtot_check=false)
    # IMPORTANT!!! ----> Mtot_target is doubled to achieve an integer representation
    channels = Int[]  # Flat vector of integers, like C++ vector<int>
    
    # Convert spin from Julia convention (actual spin value) to C++ convention (2*spin as integer)
    # e.g., spin=0.5 -> S1=1, spin=1.0 -> S1=2, spin=6.5 -> S1=13
    S1 = S_input
    S2 = S1 
    S3 = S1
    S4 = S1

    validate_Mtot_possibility(Mtot_target, JMAX, S_input, how_many_molecules)
    
    if how_many_molecules == 1
        for j1 in 0:JMAX
            for m1 in -j1:j1
                for ms1 in -S1:2:S1  # Step by 2, just like C++: ms1 = ms1 + 2
                    if check_mtot_constraint([m1], [ms1], Mtot_target; mode=constraint_check_mode, override=override_mtot_check)
                        push!(channels, j1)
                        push!(channels, m1)
                        push!(channels, S1)
                        push!(channels, ms1)
                    end
                end
            end
        end
        
    elseif how_many_molecules == 2
        for j1 in 0:JMAX
            for j2 in 0:JMAX
                for m1 in -j1:j1
                    for m2 in -j2:j2
                        for ms1 in -S1:2:S1
                            for ms2 in -S2:2:S2
                            
                                if check_mtot_constraint([m1, m2], [ms1, ms2], Mtot_target; mode=constraint_check_mode, override=override_mtot_check)
                                    push!(channels, j1)
                                    push!(channels, m1) 
                                    push!(channels, S1)
                                    push!(channels, ms1)
                                    push!(channels, j2)
                                    push!(channels, m2)
                                    push!(channels, S2) 
                                    push!(channels, ms2)
                                end
                            end
                        end
                    end
                end
            end
        end
        
    elseif how_many_molecules == 3
        for j1 in 0:JMAX
            for j2 in 0:JMAX
                for j3 in 0:JMAX
                    for m1 in -j1:j1
                        for m2 in -j2:j2
                            for m3 in -j3:j3
                                for ms1 in -S1:2:S1
                                    for ms2 in -S2:2:S2
                                        for ms3 in -S3:2:S3
                                            if check_mtot_constraint([m1, m2, m3], [ms1, ms2, ms3], Mtot_target; mode=constraint_check_mode, override=override_mtot_check)
                                                push!(channels, j1)
                                                push!(channels, m1)
                                                push!(channels, S1)
                                                push!(channels, ms1)
                                                push!(channels, j2)
                                                push!(channels, m2)
                                                push!(channels, S2)
                                                push!(channels, ms2)
                                                push!(channels, j3)
                                                push!(channels, m3)
                                                push!(channels, S3)
                                                push!(channels, ms3)
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
       
    elseif how_many_molecules == 4
        for j1 in 0:JMAX
            for j2 in 0:JMAX
                for j3 in 0:JMAX
                    for j4 in 0:JMAX
                        for m1 in -j1:j1
                            for m2 in -j2:j2
                                for m3 in -j3:j3
                                    for m4 in -j4:j4
                                        for ms1 in -S1:2:S1
                                            for ms2 in -S2:2:S2
                                                for ms3 in -S3:2:S3
                                                    for ms4 in -S4:2:S4
                                                        if check_mtot_constraint([m1, m2, m3, m4], [ms1, ms2, ms3, ms4], Mtot_target; mode=constraint_check_mode, override=override_mtot_check)
                                                            push!(channels, j1)
                                                            push!(channels, m1)
                                                            push!(channels, S1)
                                                            push!(channels, ms1)
                                                            push!(channels, j2)
                                                            push!(channels, m2)
                                                            push!(channels, S2)
                                                            push!(channels, ms2)
                                                            push!(channels, j3)
                                                            push!(channels, m3)
                                                            push!(channels, S3)
                                                            push!(channels, ms3)
                                                            push!(channels, j4)
                                                            push!(channels, m4)
                                                            push!(channels, S4)
                                                            push!(channels, ms4)
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    QN = how_many_molecules * QN_per_mol
    basis_size = length(channels) ÷ QN
    
    if print_debug
        println("Basis size: $basis_size")
        
        # Print header like C++
        if how_many_molecules == 1
            println("it\tj1\tm1\tS1\tms1")
        elseif how_many_molecules == 2  
            println("it\tj1\tm1\tS1\tms1\tj2\tm2\tS2\tms2")
        elseif how_many_molecules == 3
            println("it\tj1\tm1\tS1\tms1\tj2\tm2\tS2\tms2\tj3\tm3\tS3\tms3")
        end
        
        # Print channels like C++ 
        channel_num = 1
        for i in 1:QN:length(channels)
            print("$channel_num\t")
            for j in 0:QN-1
                print("$(channels[i+j])\t")
            end
            println()
            channel_num += 1
        end
        println()
    end
    
    return channels
end

"""
    cpp_channels_to_julia_format(cpp_channels::Vector{Int}, how_many_molecules::Int)

Convert C++ style flat integer vector to Julia style vector of Float64 vectors.
This allows comparison between the two implementations.

# Arguments
- `cpp_channels::Vector{Int}`: Flat vector of quantum numbers from C++ faithful implementation
- `how_many_molecules::Int`: Number of molecules to determine channel structure

# Returns
- `Vector{Vector{Float64}}`: Julia format channels as vector of [j,m,S,ms] vectors per molecule
"""
function cpp_channels_to_julia_format(cpp_channels::Vector{Int}, how_many_molecules::Int, QN_per_mol::Int=4)
    QN = how_many_molecules * QN_per_mol
    basis_size = length(cpp_channels) ÷ QN
    
    julia_channels = Vector{Vector{Float64}}()
    
    for i in 1:basis_size
        channel = Float64[]
        base_idx = (i-1) * QN
        
        for mol in 1:how_many_molecules
            mol_base = base_idx + (mol-1) * QN_per_mol
            j = Float64(cpp_channels[mol_base + 1])
            m = Float64(cpp_channels[mol_base + 2]) 
            S_cpp = cpp_channels[mol_base + 3]
            ms_cpp = cpp_channels[mol_base + 4]
            
            # Convert back to Julia conventions
            S_julia = Float64(S_cpp) / 2.0  # C++ stores 2*spin as integer
            ms_julia = Float64(ms_cpp) / 2.0  # C++ stores 2*ms as integer
            
            push!(channel, j)
            push!(channel, m)
            push!(channel, S_julia)
            push!(channel, ms_julia)
        end
        
        push!(julia_channels, channel)
    end
    
    return julia_channels
end

"""
    generate_channels_cpp_style(Mtot::Float64, print::Bool=false; JMAX::Int=1, spin::Float64=0.5, geometry::Vector{Float64}=Float64[])

Convert Julia format to C++ compatible parameters and call the faithful port.
This provides a Julia-compatible interface to the C++ faithful implementation.

# Arguments
- `Mtot::Float64`: Target total magnetic quantum number (Julia convention)
- `print::Bool=false`: Whether to print debug output
- `JMAX::Int=1`: Maximum rotational quantum number
- `spin::Float64=0.5`: Spin quantum number (Julia convention)
- `geometry::Vector{Float64}=Float64[]`: Molecular geometry for multi-molecule systems

# Returns
- `Vector{Vector{Float64}}`: Julia format channels compatible with generate_channels()
"""
function generate_channels_cpp_style(Mtot::Float64, print_debug::Bool=false; JMAX::Int=1, spin::Float64=0.5, geometry::Vector{Float64}=Float64[], QN_per_mol::Int=4, constraint_check_mode="simple sum", generate_all=false)
    # Convert to C++ style parameters
    Mtot_doubled = round(Int, 2 * Mtot)  # Convert Mtot to integer (doubling needed for integer operations in basis generation)
    S_cpp = round(Int, 2 * spin)  # Convert spin to C++ convention (2*S)
    how_many_molecules = isempty(geometry) ? 1 : length(geometry) ÷ 2 + 1  # Consistent calculation
    
    # Call faithful C++ port
    cpp_channels = generate_channels_cpp_faithful(Mtot_doubled, print_debug; JMAX=JMAX, S_input=S_cpp, how_many_molecules=how_many_molecules, QN_per_mol=QN_per_mol, constraint_check_mode=constraint_check_mode, override_mtot_check=generate_all)

    # Convert back to Julia format
    return cpp_channels_to_julia_format(cpp_channels, how_many_molecules, QN_per_mol)
end

end # module BasisCppFaithful
