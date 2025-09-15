# Global parameters and flags will go here.
module Parameters

export init_from_args, PARAMS
export D, Debye, ElectricDipole # Export the unit and dimension for use in other modules

using Unitful
using UnitfulAtomic
using Unitful: °, kV, cm, G, kHz, MHz, Hz, μB, GHz, nm, rad, V, s, m, kg, J, K, A, mol, Quantity # Explicitly import Quantity

# Define the electric dipole dimension if not already defined
@derived_dimension ElectricDipole Unitful.𝐋*Unitful.𝐈*Unitful.𝐓

# Define the Debye unit (D) in terms of elementary charge and Bohr radius
# 1 D ≈ 0.393456 e⋅a₀ in atomic units as per the C++ implementation
# Note: e_au and a0_au are the elementary charge and Bohr radius in atomic units
@unit D "D" Debye 0.393456*UnitfulAtomic.e_au*UnitfulAtomic.a0_au false

# Register the unit to make it available to the u"" string macro
Unitful.register(Parameters)

# Make sure the Debye unit is preferred for electric dipole dimensions
if !@isdefined(ElectricDipoleUnits)
    const ElectricDipoleUnits = Union{typeof(u"D")}
    Unitful.preferunits(u"D")
    
    # Add promotion rule for ElectricDipole dimension to prefer Debye
    Unitful.promote_unit(::S, ::T) where {S<:ElectricDipoleUnits, T<:ElectricDipoleUnits} = u"D"
end


# Ensure specific units are available for constructing quantities if not fully qualified.
# Alternatively, use u"string" macro for all.
# using Unitful: kV, cm, °, G, kHz, MHz, Hz, D, μB, GHz, nm # Already present

mutable struct Params
    # Core calculation mode
    mode::Int # 0 = energy, 1 = magnetization, 2 = time evolution, etc. (changed to Float64 to match parse usage)
    vsWhat::Int # This might determine which parameter is being varied in a scan
    
    # Geometry - Stores (R, Phi) for each additional molecule relative to first one.
    # R in nm, Phi in degrees from vertical below first molecule (0 = below, 90 = right, etc)
    # The first molecule is always at (0,0). Each argument adds a molecule relative to the first.
    # Only 2 or 3 molecules supported (so 2 or 4 values: (R1, Phi1), (R2, Phi2))
    interaction_geometries::Vector{Tuple{Quantity{Float64, dimension(u"nm")}, Quantity{Float64, dimension(u"°")}}}
    how_many_molecules::Int  # Total number of molecules (calculated from geometry)

    # External fields - Store as Unitful quantities
    E::Quantity{Float64, dimension(u"kV/cm")}
    Phi_el::Quantity{Float64, dimension(u"°")}
    Theta_el::Quantity{Float64, dimension(u"°")}
    
    B::Quantity{Float64, dimension(u"Gauss")}
    Phi_mg::Quantity{Float64, dimension(u"°")}
    Theta_mg::Quantity{Float64, dimension(u"°")}
    
    # Molecular properties - Store as Unitful quantities
    A::Quantity{Float64, dimension(u"kHz")} 
    d_el::Quantity{Float64, dimension(u"D")} 
    d_mg::Quantity{Float64, dimension(u"μB")}
    
    # Iteration parameters
    Mtot::Float64 
    START_IT::Int
    END_IT::Int 
    # step can be unitful if it represents a physical quantity being scanned
    # For now, assume it's a Float64, but could be made more generic or Unitful.
    step::Float64
    IT_NO_OVERRIDE::Union{Nothing, Int}
    ITS_PER_OVERRIDE::Union{Nothing, Int}

    # Basis parameters
    JMAX::Int
    QN_per_mol::Int  # Number of quantum numbers per molecule (default: 4 for j,m,S,ms)
    Brot::Quantity{Float64, dimension(u"GHz")}
    spin::Float64 # Input is 2S
                               
    # Quench parameters
    quenched_E::Union{Nothing, Quantity{Float64, dimension(u"kV/cm")}}
    quenched_B::Union{Nothing, Quantity{Float64, dimension(u"Gauss")}}
    quenched_R::Union{Nothing, Vector{Tuple{typeof(1.0u"nm"), typeof(1.0u"°")}}}
    observables::Vector{String}

    # Product state evolution
    productstate::Union{Nothing, String}
    initial_eigenstate::Union{Nothing, Int}
    initial_basis_state::Union{Nothing, Int}
    initial_state_label::Union{Nothing, String}  # Label for the initial state, if applicable
    
    # Harmonic time parameters
    harmonic_frequency::Union{Nothing, Quantity{Float64}}  # The ω for harmonic units
    time_is_harmonic::Bool  # Flag to indicate time is in harmonic units
    tmax::Float64    # Time in harmonic units (dimensionless)
    dt::Float64      # Time step in harmonic units (dimensionless)

    # Output directory
    outputdir::Union{Nothing, String}
    
    # Flags
    MG_QUENCH::Bool

    # Data type for Hamiltonian calculations
    hamiltonian_dtype::String # e.g., :Float64, :ComplexF64, etc.
    spectra_unit::String # In what units are spectra saved 

    # Skip terms
    skip_terms::Union{Nothing, Set{Symbol}}

    geom_values_flat::Vector{Float64}

    # Looped eigenproblem flag
    eigenproblem_looped::Bool

    # Flag for colored printing in terminal
    print_colored::Bool
    progressbar::Bool

    # Run_id_string
    run_id_string::Union{Nothing, String}

    # Profiling enabled
    profile::Bool

    # Root folder for saving results
    save_root_folder::String # Default is current directory, can be set to scratch on cluster


    # Constructor with default Unitful quantities
    Params() = new(
        0.0, 0, # mode, vsWhat
        Vector{Tuple{Quantity{Float64, dimension(u"nm")}, Quantity{Float64, dimension(u"°")}}}(), # Initialize interaction_geometries as empty
        1, # how_many_molecules default (1 molecule if no geometry)
        0.0u"kV/cm", 0.0u"°", 0.0u"°",
        0.0u"Gauss", 0.0u"°", 0.0u"°",
        0.0u"kHz", 0.0u"D", 0.0u"μB",
        0.0, 0, 0, 1.0, nothing, nothing, # Mtot, START_IT, END_IT, step, IT_NO_OVERRIDE, ITS_PER_OVERRIDE
        0, 4, 0.0u"GHz", 0.0, # JMAX, QN_per_mol, Brot, spin
        nothing, nothing, nothing, # quenched_E, quenched_B, quenched_R
        Vector{String}(), # observables
        nothing, 0, nothing, nothing, # productstate, initial_basis_state, initial_eigenstate, initial_state_label
        nothing, false, 1.0, 0.01,  # harmonic_frequency, time_is_harmonic, tmax, dt
        nothing, # outputdir
        false,    # MG_QUENCH
        "ComplexF64", # hamiltonian_dtype
        "atomic", # Spectra unit
        nothing, # skip_terms
        Vector{Float64}(), # geom_values_flat
        false, # If eigenproblem is looped
        false, # Colored printing flag
        false, # Progress bar off by default
        nothing, # In default do not add parallel run ids
        false, # Profiling off by default
        "./", # Default save root folder is current directory
    )
end

const PARAMS = Ref{Union{Nothing, Params}}(nothing)

# Helper to safely parse units, returning a default or erroring
function parse_unitful_arg(value::Real, unit_str::String, default_unit::Unitful.Units)
    parsed_unit = try
        uparse(unit_str)
    catch e
        @warn "Could not parse unit string '$unit_str'. Using default unit '$default_unit'. Error: $e"
        default_unit # Or rethrow, or handle more gracefully
    end
    if dimension(parsed_unit) != dimension(default_unit)
        error("Parsed unit '$parsed_unit' for value '$value' is not compatible with expected dimension of '$default_unit'.")
    end
    return value * parsed_unit
end

function get_skip_term_for_vsWhat(; vsWhat::Int=-1, str_skip_term="")
    if vsWhat == 0 || str_skip_term == "electric_field"
        return [:electric_field]
    elseif vsWhat == 1 || str_skip_term == "magnetic_field"
        return [:magnetic_field]
    elseif vsWhat == 2 || str_skip_term == "spin_rotation"
        return [:spin_rotation]
    elseif vsWhat == 3 || str_skip_term == "geometry_scaling"
        return [:electric_dipolar, :magnetic_dipolar]  # For geometry scaling, we rebuild dipolar interactions
    elseif vsWhat == 4 || str_skip_term == "electric_dipolar"
        return [:electric_dipolar]  # When iterating d_el, skip dipolar interactions in base
    elseif vsWhat == 5 || str_skip_term == "magnetic_dipolar"
        return [:magnetic_dipolar]  # When iterating d_mg, skip magnetic dipolar in base
    else
        error("Unknown vsWhat: $vsWhat")
    end
end

function init_from_args(args)
    # Remove all config file merging logic
    current_params = Params()
    mode_vsWhat = string(args["mode"])
    if occursin('.', mode_vsWhat)
        # Parse as float first, then extract integer and fractional parts
        mode_float = parse(Float64, mode_vsWhat)
        current_params.mode = Int(floor(mode_float))  # Integer part
        fractional_part = mode_float - floor(mode_float)
        current_params.vsWhat = Int(round(fractional_part * 10))  # Fractional part * 10
    else
        current_params.mode = parse(Int, mode_vsWhat)
        current_params.vsWhat = 0
    end
    current_params.hamiltonian_dtype = args["hamiltonian-dtype"] # Assuming this is a string or symbol

    # Parse geometry: expect 2 or 4 values (for 2 or 3 molecules)
    geom_values_flat = args["geometry..."]::Vector{Float64}
    current_params.geom_values_flat = geom_values_flat
    if !(length(geom_values_flat) == 0 || length(geom_values_flat) == 2 || length(geom_values_flat) == 4)
        error("Geometry must be 0 (single molecule), 2 (two molecules), or 4 (three molecules): (R1, Phi1), (R2, Phi2). Got $(length(geom_values_flat)) values.")
    end
    current_params.interaction_geometries = Tuple{Quantity{Float64, dimension(u"nm")}, Quantity{Float64, dimension(u"°")}}[]
    for i in 1:2:length(geom_values_flat)
        R = geom_values_flat[i] * u"nm"
        Phi = geom_values_flat[i+1] * u"°"
        push!(current_params.interaction_geometries, (R, Phi))
    end

    current_params.how_many_molecules = isempty(geom_values_flat) ? 1 : length(geom_values_flat) ÷ 2 + 1

    # Electric Field
    E_val = args["E"]
    E_unit_str = args["E-unit"]
    current_params.E = parse_unitful_arg(E_val, E_unit_str, u"kV/cm")
    current_params.Phi_el = args["phi_el"] * u"°" # CLI gives degrees directly
    current_params.Theta_el = args["theta_el"] * u"°"

    # Magnetic Field
    current_params.B = args["B"] * u"Gauss"
    current_params.Phi_mg = args["phi_mg"] * u"°"
    current_params.Theta_mg = args["theta_mg"] * u"°"

    # Spin-Rotation Coupling
    A_val = args["A"]
    A_unit_str = args["A-unit"]
    current_params.A = parse_unitful_arg(A_val, A_unit_str, u"kHz")

    # Dipole Moments
    current_params.d_el = args["d_el"] * u"D" # Assuming CLI gives Debye directly
    current_params.d_mg = args["d_mg"] * u"μB" # Assuming CLI gives μB directly

    # Dimensionless/Counts
    current_params.Mtot = args["Mtot"]
    current_params.START_IT = args["START_IT"]
    current_params.END_IT = args["END_IT"]
    current_params.step = args["step"] # Assuming dimensionless or to be handled later if unitful
    if args["IT_NO_OVERRIDE"] !== nothing
        current_params.IT_NO_OVERRIDE = args["IT_NO_OVERRIDE"]
    end

    if args["ITS_PER_OVERRIDE"] !== nothing
        current_params.ITS_PER_OVERRIDE = args["ITS_PER_OVERRIDE"]
    end

    current_params.JMAX = args["JMAX"]
    current_params.QN_per_mol = get(args, "QN_per_mol", 4)  # Default to 4 if not provided
    current_params.spin = args["spin"] # This is 2S

    # Rotational Constant
    Brot_val = args["Brot"]
    Brot_unit_str = args["Brot-unit"]
    current_params.Brot = parse_unitful_arg(Brot_val, Brot_unit_str, u"GHz")

    # Quenched Parameters (handle potential nothing)
    if args["quenched_E"] !== nothing
        current_params.quenched_E = args["quenched_E"] * u"kV/cm" 
    end
    if args["quenched_B"] !== nothing
        current_params.quenched_B = args["quenched_B"] * u"Gauss"
    end
    if args["quenched_R"] !== nothing
        new_values_flat = args["quenched_R"]::Vector{Float64}
        new_interaction_geometries = Tuple{Quantity{Float64, dimension(u"nm")}, Quantity{Float64, dimension(u"°")}}[]

        for i in 1:2:length(new_values_flat)
            R = new_values_flat[i] * u"nm"
            Phi = new_values_flat[i+1] * u"°"
            push!(new_interaction_geometries, (R, Phi))
        end
        current_params.quenched_R = new_interaction_geometries
    end
    
    # Determine harmonic frequency
    harmonic_freq_quantity = if args["harmonic-frequency-source"] == "Brot"
        current_params.Brot
    elseif args["harmonic-frequency-source"] == "A"
        current_params.A
    elseif args["harmonic-frequency-source"] == "custom"
        if args["harmonic-frequency"] === nothing
            error("Custom harmonic frequency source selected but --harmonic-frequency not provided")
        end
        freq_val = args["harmonic-frequency"]
        freq_unit_str = args["harmonic-frequency-unit"]
        parse_unitful_arg(freq_val, freq_unit_str, u"MHz") * (2π) # Added (2π) to convert to angular frequency
    else
        error("Invalid harmonic frequency source: $(args["harmonic-frequency-source"])")
    end

    # Store harmonic frequency in parameters
    current_params.harmonic_frequency = harmonic_freq_quantity

    current_params.observables = get(args, "quench-observables", ["Ms", "M"])

    # Handle time parameters based on time unit
    time_unit = args["time-unit"]
    tmax_val = args["tmax"]
    dt_val = args["dt"]

    if time_unit == "harmonic"
        # Input is already in harmonic units (dimensionless)
        current_params.tmax = tmax_val  # dimensionless harmonic time
        current_params.dt = dt_val      # dimensionless harmonic time
        current_params.time_is_harmonic = true
    else
        # Input is in SI units, need to convert to harmonic
        try
            tmax_quantity = tmax_val * uparse(time_unit)
            dt_quantity = dt_val * uparse(time_unit)
            
            # Convert SI time to harmonic units using the harmonic frequency
            current_params.tmax = time_si_to_harmonic_units(tmax_quantity, harmonic_freq_quantity)
            current_params.dt = time_si_to_harmonic_units(dt_quantity, harmonic_freq_quantity)
            current_params.time_is_harmonic = true
            
            println("Converted time from SI to harmonic units:")
            println("  tmax: $(tmax_val) $time_unit → $(current_params.tmax) × (2π/ω)")
            println("  dt: $(dt_val) $time_unit → $(current_params.dt) × (2π/ω)")
            
        catch e
            error("Invalid time unit '$time_unit'. Supported: s, ms, μs, ns, ps, fs, or 'harmonic'. Error: $e")
        end
    end
    
    # Initial state
    current_params.productstate = get(args, "productstate", nothing)

    # Initial eigenstate
    # Assuming initial_eigenstate is an integer, default to 1 if not provided
    current_params.initial_eigenstate = get(args, "initial-eigenstate", nothing)

    # Assuming initial_basis_state is an integer, default to 1 if not provided
    current_params.initial_basis_state = get(args, "initial-basis-state", nothing)

    current_params.initial_state_label = nothing

    current_params.outputdir = args["outputdir"]

    current_params.eigenproblem_looped = args["looped-eigenproblem"]

    current_params.spectra_unit = args["spectra-unit"]

    skip_terms_str = get(args, "skip-terms", nothing)
    if skip_terms_str !== nothing
        skip_terms_clean = [strip(term) for term in split(skip_terms_str, ',') if !isempty(strip(term))]
        skip_terms_symbols = [get_skip_term_for_vsWhat(str_skip_term=term) for term in skip_terms_clean]
        skip_terms_symbols = [(skip_terms_symbols...)...]
        println("Artificial skip terms detected: $skip_terms_symbols")
        current_params.skip_terms = Set(skip_terms_symbols)
    else
        current_params.skip_terms = nothing
    end

    current_params.print_colored = get(args, "colored-print", false)

    current_params.progressbar = args["progressbar"]

    if args["run-id-string"] !== nothing
        current_params.run_id_string = args["run-id-string"]
    end

    current_params.save_root_folder = get(args, "save-root-folder", "./") # Default to current directory

    PARAMS[] = current_params
end

end # module Parameters
