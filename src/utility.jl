# Utility functions: file naming, saving, printing, etc.
module Utility

export create_filename, get_quench_filename, get_data_folder, params_to_cpp_strings

using ..Parameters
using Unitful: ustrip
using Unitful
using LinearAlgebra
using SparseArrays
using Arrow
using JLD2

# Helper function to parse dtype string to actual type
function parse_hamiltonian_dtype(dtype_str::String)
    if dtype_str == "ComplexF64"
        return ComplexF64
    elseif dtype_str == "ComplexF32"
        return ComplexF32
    elseif dtype_str == "Float64"
        return Float64
    elseif dtype_str == "Float32"
        return Float32
    elseif dtype_str == "ComplexF16"
        return ComplexF16
    else
        error("Unsupported Hamiltonian data type: $dtype_str")
    end
end

function get_loop_value(params, iterator)
    # For each vsWhat, return the value to use for the iterated parameter at this iterator
    # iterator==0 should correspond to the value in params
    val = if params.vsWhat == 0
        # Electric field
        ustrip(params.E) + iterator * params.step
    elseif params.vsWhat == 1
        # Magnetic field
        ustrip(params.B) + iterator * params.step
    elseif params.vsWhat == 2
        # Spin-rotation
        ustrip(params.A) + iterator * params.step
    elseif params.vsWhat == 3
        # Geometry scaling
        1.0 + iterator * params.step
    elseif params.vsWhat == 4
        # Electric dipole
        ustrip(params.d_el) + iterator * params.step
    elseif params.vsWhat == 5
        # Magnetic dipole
        ustrip(params.d_mg) + iterator * params.step
    else
        error("Unknown vsWhat: $(params.vsWhat)")
    end
    return round(val, digits=5)
end

function get_data_folder(mode::Int, subtype::Int=0;params_override::Union{Nothing, Parameters.Params}=nothing)
    intpart = floor(Int, mode)
    
    if params_override !== nothing
        # Use provided params override if available
        params = params_override
    else
        # Get parameters from Parameters module
        params = Parameters.PARAMS[]
    end
    if params === nothing
        error("Parameters not initialized when trying to get data folder.") # Added error handling
    end

    # NEW: Add molecule count folder
    how_many_molecules = getfield(params, :how_many_molecules)
    mol_label = how_many_molecules == 1 ? "1_molecule" : "$(how_many_molecules)_molecules"
    
    # Extract parameter values for directory hierarchy
    # Use ustrip to get the numerical value from the Unitful.Quantity
    E_field_val = ustrip(params.E)
    B_field_val = ustrip(params.B)
    A_factor_val = ustrip(params.A)
    D_EL_val = params.d_el / u"D"  # Extract value in Debye
    D_MG_val = params.d_mg / u"μB"  # Extract value in Bohr magnetons
    
    # Create 5-level parameter-based subdirectory structure like C++
    # Using Int for integer-like values, otherwise format as needed.
    # Consider a helper function for consistent string formatting of these values if more control is needed.
    EF_dir = "EF$(Int(round(E_field_val)))" # Rounding before Int in case of small float inaccuracies if inputs were float
    MF_dir = "MF$(Int(round(B_field_val)))"
    A_dir = "A$(Int(round(A_factor_val)))"
    # For D_EL, which might be float, decide on formatting. Using a simple approach for now.
    D_EL_dir = "D_EL$(D_EL_val == floor(D_EL_val) ? Int(D_EL_val) : D_EL_val)" 
    D_MG_dir = "D_MG$(Int(round(D_MG_val)))"
    
    base_data_path = "data" # Base directory for all data

    # I think this accidentally ommits the EF level in folder structure determination in quench saving (vsWhat=0)
    if params.mode in [0, 2]
        if params.skip_terms === nothing
            skip_terms = Parameters.get_skip_term_for_vsWhat(vsWhat=params.vsWhat)
        else
            skip_terms_tmp = Parameters.get_skip_term_for_vsWhat(vsWhat=params.vsWhat)
            skip_terms = union(skip_terms_tmp, params.skip_terms)
        end
    else
        skip_terms = params.skip_terms === nothing ? Set{Symbol}() : params.skip_terms
    end

    levels_list = []

    for (name, level) in [(:electric_field, EF_dir), (:magnetic_field, MF_dir), (:spin_rotation, A_dir), (:electric_dipolar, D_EL_dir), (:magnetic_dipolar, D_MG_dir)]
        if name in skip_terms
            # println("Iteration over: $name")
            continue
        end
        push!(levels_list, level)
    end
    if params.eigenproblem_looped
        iterator_label = get_vs_folder_name(params.vsWhat)[3:end]
        val_min = format_number_string(get_loop_value(params,0))
        val_max = format_number_string(get_loop_value(params,params.END_IT-1))
        run_id = "$(iterator_label)min$(val_min)_$(iterator_label)max$(val_max)_d$(iterator_label)$(params.step)"
        if params.run_id_string !== nothing
            run_id *= "_" * params.run_id_string
        end
        push!(levels_list, run_id)
    end

    if intpart == 0 || intpart == 2 # Spectrum calculations
        # Define path for spectrum results, e.g.,
        return joinpath(params.save_root_folder, base_data_path, "spectrum_results", mol_label, levels_list...)
    elseif intpart == 1  # Quench dynamics
        return joinpath(params.save_root_folder, base_data_path, "quench_results", mol_label, levels_list...)
    # elseif intpart == 2 # Product state evolution
    #     return joinpath(params.save_root_folder, base_data_path, "product_state_results", mol_label, levels_list...)
    elseif intpart == 3 # Magnetization
        return joinpath(params.save_root_folder, base_data_path, "magnetization_results", mol_label, levels_list...)
    elseif intpart == 4 # Eigenfunctions
        return joinpath(params.save_root_folder, base_data_path, "eigenfunction_results", mol_label, levels_list...)
    else
        error("Unknown mode $mode for data folder generation.")
    end
end

function create_filename(
    params::Parameters.Params, 
    iterator, 
    what_quantity::String; 
    override_top_folder::Union{Nothing,String}=nothing, 
    file_extension::String=".txt",
    override_params=false,
    )
    # Geometry part
    geometry_length = 2 * (params.how_many_molecules - 1)
    string_geometry = "geo"
    for i in 1:geometry_length
        string_geometry *= "_" * format_number_string.(params.geom_values_flat[i])
    end
    string_iterator = lpad(string(iterator), 3, '0')
    filename_components = [
        string_iterator,
        what_quantity,
        "Jmax$(  format_number_string(ustrip(params.JMAX)))",
        "Brot$(  format_number_string(ustrip(params.Brot)))",
        "s$(     format_number_string(ustrip(params.spin)))",
        "E$(     format_number_string(ustrip(params.E)))",
        "thetaE$(format_number_string(ustrip(params.Theta_el)))",
        "phiE$(  format_number_string(ustrip(params.Phi_el)))",
        "B$(     format_number_string(ustrip(params.B)))",
        "thetaB$(format_number_string(ustrip(params.Theta_mg)))",
        "phiB$(  format_number_string(ustrip(params.Phi_mg)))",
        "A$(     format_number_string(ustrip(params.A)))",
        "del$(   format_number_string(ustrip(params.d_el / u"D")))",
        "dmg$(   format_number_string(ustrip(params.d_mg / u"μB")))",
        "Mtot$(  format_number_string(ustrip(params.Mtot)))",
        string_geometry*file_extension
    ]
    filename = join(filename_components, "_")

    # Get appropriate data folder and ensure it exists
    mode = params.mode
    subtype = floor(Int, (mode - floor(mode)) * 10)
    if override_params
        data_folder = get_data_folder(mode, subtype, params_override=params)
    else
        data_folder = get_data_folder(mode, subtype)
    end


    # If override_top_folder is given, replace the topmost folder *within* spectrum_results with it
    if override_top_folder !== nothing
        parts = splitpath(data_folder)
        idx = findfirst(x -> x == "spectrum_results", parts)
        if idx !== nothing && idx < length(parts)
            # Replace the next part (the folder within spectrum_results) with override_top_folder
            parts[idx+1] = override_top_folder
            data_folder = joinpath(parts...)
        end
    end
    
    mkpath(data_folder)

    return joinpath(data_folder, filename)
end

function get_quench_filename(
    current_params::Parameters.Params,
    observable::String,
    save_folder::String;
    couplings=false,
    consolidated=false,
    )
    if couplings
        filename = joinpath(save_folder, "quench_$(observable)_couplings.csv")
        return filename
    end
    if current_params.IT_NO_OVERRIDE !== nothing && !consolidated
        if current_params.ITS_PER_OVERRIDE !== nothing
            filename = joinpath(save_folder, "IT$(current_params.IT_NO_OVERRIDE)-$(current_params.IT_NO_OVERRIDE+current_params.ITS_PER_OVERRIDE - 1)_quench_$(observable).csv")
        else
            filename = joinpath(save_folder, "IT$(current_params.IT_NO_OVERRIDE)_quench_$(observable).csv")
        end
    else
        filename = joinpath(save_folder, "quench_observable_$(observable).csv")
    end
    return filename
end

# Helper function to format numbers as integers when they have no decimal part
function format_number_string(value::Real)
    if value == floor(value)
        return string(Int(value))
    else
        return string(value)
    end
end

# Convert parameters to C++-compatible string format for filename generation
# This function strips units and formats parameters exactly as the C++ version expects
function params_to_cpp_strings(params)
    # Create the parameter strings vector in the exact order expected by create_filename
    # Order: geometry, E, thetaE, phiE, B, thetaB, phiB, A, d_el, d_mg, Mtot, JMAX, B_rot, spin, Rx
    
    string_params = String[]
    
    # Add geometry parameters (already in correct format)
    for geo_param in params.geom_values_flat
        push!(string_params, format_number_string(geo_param))
    end
    
    # Electric field parameters
    push!(string_params, format_number_string(ustrip(params.E)))  # Strip kV/cm unit
    push!(string_params, format_number_string(ustrip(params.Theta_el)))   # Already in degrees
    push!(string_params, format_number_string(ustrip(params.Phi_el)))     # Already in degrees

    # Magnetic field parameters
    push!(string_params, format_number_string(ustrip(params.B)))  # Strip Gauss unit
    push!(string_params, format_number_string(ustrip(params.Theta_mg)))   # Already in degrees
    push!(string_params, format_number_string(ustrip(params.Phi_mg)))     # Already in degrees

    # Spin-rotation coupling
    push!(string_params, format_number_string(ustrip(params.A)))  # Strip kHz unit
    
    # Electric dipole moment (convert to Debye)
    d_el_debye = ustrip(params.d_el / u"D")
    push!(string_params, format_number_string(d_el_debye))
    
    # Magnetic dipole moment (convert to Bohr magnetons)
    d_mg_bohr = ustrip(params.d_mg / u"μB")
    push!(string_params, format_number_string(d_mg_bohr))
    
    # Total magnetic quantum number
    push!(string_params, format_number_string(params.Mtot))
    
    # Basis parameters
    push!(string_params, format_number_string(params.JMAX))
    push!(string_params, format_number_string(ustrip(params.Brot)))  # Strip GHz unit
    push!(string_params, format_number_string(params.spin))

    
    return string_params
end

# Get the folder for saving quench data, based on mode, subtype, and quench-specific labels
function get_quench_save_folder(current_params::Parameters.Params; override=false)
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

    time_unit = "harmonic"  # Use CLI time unit directly
    time_suffix = "tmax_$(current_params.tmax)dt_$(current_params.dt)unit_$(time_unit)"
    
    if override 
        base_folder = get_data_folder(current_params.mode, current_params.vsWhat, params_override=current_params)
    else
        base_folder = get_data_folder(current_params.mode, current_params.vsWhat)
    end
    return joinpath(base_folder, quench_suffix, current_params.initial_state_label, time_suffix)
end

function get_vs_folder_name(vsWhat::Int)
    if vsWhat == 0
        return "vsE"
    elseif vsWhat == 1
        return "vsB"
    elseif vsWhat == 2
        return "vsA"
    elseif vsWhat == 3
        return "vsR"
    elseif vsWhat == 4
        return "vsDEL"
    elseif vsWhat == 5
        return "vsDMG"
    else
        return "LOOPED"
    end
end

function get_var_dtype(dtype_string)
    if dtype_string == "Float64" || dtype_string == "ComplexF64"
        return Float64
    elseif dtype_string == "Float32" || dtype_string == "ComplexF32"
        return Float32
    elseif dtype_string == "Float16"
        return Float16
    else
        error("Wrong dtype")
    end
end


function get_spectrum_filename(
    params::Parameters.Params;
    iterator="",
    object_name="object",
    subfolder::String="",
    just_dir::Bool=false,
    file_extension=".arrow",
    override_params=false,
)
    mol_label = params.how_many_molecules == 1 ? "1_molecule" : "$(params.how_many_molecules)_molecules"
    override_top_folder = joinpath(mol_label, get_vs_folder_name(params.vsWhat))
    filename = create_filename(params, iterator, object_name; override_top_folder=override_top_folder, file_extension=file_extension, override_params=override_params)
    if !isempty(subfolder)
        folder = dirname(filename)
        filename_only = basename(filename)
        folder_with_sub = joinpath(folder, subfolder)
        mkpath(folder_with_sub)
        filename = joinpath(folder_with_sub, filename_only)
    else
        dir_path = dirname(filename)
        if !isdir(dir_path)
            mkpath(dir_path)
        end
    end

    if just_dir
        return dirname(filename)
    end
    return filename
end

function get_filename(
    params::Parameters.Params;
    object_name="object",
    override_params=false,
    just_dir=false,
    ensure_mkdir=true,
    get_filename=false,
    kwargs...
)
    if params.mode in [0, 2]
        # Spectrum calculation mode
        filename = get_spectrum_filename(params; object_name=object_name, override_params=override_params, just_dir=just_dir, kwargs...)
        if !just_dir
            folder = dirname(filename)
            if ensure_mkdir && !isdir(folder)
                mkpath(folder)
            end
        else
            if ensure_mkdir && !isdir(filename)
                mkpath(filename)
            end
        end
        return filename
    elseif params.mode == 1
        if just_dir
            return get_quench_save_folder(params; override=override_params)
        end
        folder = get_quench_save_folder(params; override=override_params)
        if ensure_mkdir && !isdir(folder)
            mkpath(folder)
        end
        filename = get_quench_filename(params, object_name, folder; kwargs...)
        return filename
    end
end

end # module
