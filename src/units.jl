# Unit conversion functions will go here.
module Units

using ..PhysicalConstantsJL # To access ħ, c_0 etc. from CODATA2018
using Unitful
using UnitfulAtomic

# Import specific units and atomic units for clarity and direct use
using Unitful: nm, cm, V, kV, MV, GV, Gauss, Hz, kHz, MHz, GHz, J, K, s, m, kg, A, mol, rad, °, Quantity, uconvert, ustrip, dimension, uparse, @u_str
using UnitfulAtomic: hartree, bohr, aunit # Explicitly import atomic unit symbols if preferred

export to_au_value, to_radians_value # Exporting functions that return values
export convert_param_to_au # Higher-level conversion for specific param types
export frequency_to_energy_au_value, wavenumber_to_energy_au_value, eigenvalues_au_to_harmonic_units # Export specialized converters

# --- Target Atomic Units Dictionary ---
# Defines the target UnitfulAtomic unit for various types of physical quantities.
const TARGET_AU_UNITS = Dict{Symbol, Unitful.Units}(
    :energy                   => u"hartree", # Hartree energy (Eh)
    :length                   => u"bohr",    # Bohr radius (a₀)
    :time                     => aunit(u"s"), # Atomic unit of time (ħ/Eh)
    :electric_dipole_moment   => UnitfulAtomic.e_au * UnitfulAtomic.a0_au, # Atomic unit of electric dipole (e*a₀)
    :magnetic_dipole_moment   => aunit(u"A*m^2"), # Atomic unit of magnetic dipole moment
    :electric_field           => aunit(u"V/m"), # Atomic unit of electric field
    :magnetic_field           => aunit(u"T"), # Atomic unit of magnetic field
    :angle                    => u"rad" # Angles are fundamentally dimensionless but often need conversion to radians
    # Add other common types as needed:
    # :mass                     => u"me_au",                 # electron mass (m_e)
    # :charge                   => u"e_au",                  # elementary charge (e)
    # :velocity                 => aunit(u"m/s"),            # bohr * hartree / ħ
    # :momentum                 => aunit(u"kg*m/s"),         # ħ / bohr
    # :angular_momentum         => u"ħ_au",                  # ħ itself
)

# --- Conversions to Atomic Unit VALUES (dimensionless Float64) ---

# Generic function to convert a Unitful quantity to its dimensionless value in a target atomic unit.
function to_au_value(quantity::Unitful.AbstractQuantity, unit_type_symbol::Symbol)
    if !haskey(TARGET_AU_UNITS, unit_type_symbol)
        error("Unknown unit type symbol for AU conversion: $unit_type_symbol. Supported: $(keys(TARGET_AU_UNITS))")
    end
    target_unit = TARGET_AU_UNITS[unit_type_symbol]
    return ustrip(uconvert(target_unit, quantity))
end

# Specific conversion wrappers/examples for convenience and clarity, returning dimensionless AU values.

# Length: e.g., nm to Bohr value
# Example: to_au_value(1.0u"nm", :length) 

# Electric Field: e.g., kV/cm to AU of electric field value
# Example: to_au_value(1.0u"kV/cm", :electric_field)

# Magnetic Field: e.g., Gauss to AU of magnetic field value
# Example: to_au_value(1.0u"G", :magnetic_field)

# Energy/Frequency: e.g., kHz, GHz, cm^-1 to Hartree value
# For frequencies (E=hf) or wavenumbers (E=hcν̃), first convert to an energy unit (e.g., Joule)
# then convert to Hartree.
function frequency_to_energy_au_value(freq_quantity::Unitful.AbstractQuantity)
    if dimension(freq_quantity) != dimension(u"Hz")
        error("Input to frequency_to_energy_au_value must be a frequency. Got: $(freq_quantity)")
    end
    # E = h * f. PhysicalConstantsJL.planck_constant is h.
    energy_joules = uconvert(u"J", freq_quantity * PhysicalConstantsJL.planck_constant)
    return to_au_value(energy_joules, :energy) # :energy maps to au_energy (Hartree)
end

function wavenumber_to_energy_au_value(wn_quantity::Unitful.AbstractQuantity)
    # Ensure wn_quantity is in a dimensionally correct unit like m^-1 or cm^-1
    if dimension(wn_quantity) != dimension(u"m^-1")
        error("Input to wavenumber_to_energy_au_value must be a wavenumber (e.g., m^-1, cm^-1). Got: $(wn_quantity)")
    end
    # E = hcν̃. PhysicalConstantsJL.planck_constant is h, PhysicalConstantsJL.speed_of_light_in_vacuum is c.
    energy_joules = uconvert(u"J", wn_quantity * PhysicalConstantsJL.planck_constant * PhysicalConstantsJL.speed_of_light_in_vacuum)
    return to_au_value(energy_joules, :energy)
end

# --- Higher-level Parameter Conversion ---

# Converts a parameter (Unitful.Quantity) from the PARAMS struct or similar
# to its dimensionless atomic unit value based on its physical nature.
function convert_param_to_au(param_value::Unitful.Quantity, param_type_symbol::Symbol)
    # param_type_symbol indicates the nature of the parameter, e.g.,
    # :Brot, :A_SR (spin-rotation), :E_field, :B_field, :d_el, :d_mg, :R_dist, :time_val
    
    if param_type_symbol === :Brot || param_type_symbol === :A_SR # Rotational constant, Spin-rotation
        # These are typically given in frequency units (e.g., GHz, kHz)
        if dimension(param_value) == dimension(u"Hz")
            return frequency_to_energy_au_value(param_value)
        elseif dimension(param_value) == dimension(u"m^-1") # If given as wavenumber (e.g., cm^-1)
            return wavenumber_to_energy_au_value(param_value)
        elseif dimension(param_value) == dimension(u"J") # If already in energy units
            return to_au_value(param_value, :energy)
        else
            error("Unsupported unit dimension for energy-like parameter '$param_type_symbol': $(dimension(param_value)). Expected frequency, wavenumber, or energy.")
        end
    elseif param_type_symbol === :E_field # Electric field strength (e.g., PARAMS[].E is kV/cm)
        return to_au_value(param_value, :electric_field)
    elseif param_type_symbol === :B_field # Magnetic field strength (e.g., PARAMS[].B is Gauss)
        return to_au_value(param_value, :magnetic_field)
    elseif param_type_symbol === :d_el # Electric dipole moment (e.g., PARAMS[].d_el is Debye)
        return to_au_value(param_value, :electric_dipole_moment)
    elseif param_type_symbol === :d_mg # Magnetic dipole moment (e.g., PARAMS[].d_mg is μB)
        # This converts the input moment (e.g., 1.0u"μB") to its value in UnitfulAtomic.au_magnetic_moment (which is eħ/m_e = 2μB).
        # So, if param_value is 1.0u"μB", this returns 0.5.
        # The g-factor (e.g., g_S=2 for electron spin) is applied separately in the Hamiltonian construction.
        return to_au_value(param_value, :magnetic_dipole_moment)
    elseif param_type_symbol === :R_dist # Inter-particle distance (e.g., geometry R is nm)
        return to_au_value(param_value, :length)
    elseif param_type_symbol === :time_val # For time parameters like tmax, dt (e.g., in seconds)
        return to_au_value(param_value, :time)
    # Add more specific parameter type conversions as the project evolves
    else
        error("Unknown parameter type for AU conversion: $param_type_symbol. Supported types depend on implementation.")
    end
end

# --- Angle Conversions ---

# Converts an angle quantity (e.g., in degrees) to its dimensionless value in radians.
function to_radians_value(angle_quantity::Unitful.AbstractQuantity)
    # Check if the quantity is an angle (dimensionless in Unitful terms, but has units like °, rad)
    if !(dimension(angle_quantity) == dimension(u"°") || dimension(angle_quantity) == dimension(u"rad") || dimension(angle_quantity) == Unitful.NoDims)
        # Allow NoDims for cases where angle might already be a bare number assumed to be in radians, though Unitful quantities are preferred.
        # Stricter check: dimension(angle_quantity) != dimension(u"rad") && dimension(angle_quantity) != dimension(u"°")
        error("Input to to_radians_value must be an angle quantity (e.g., in degrees or radians). Got: $(angle_quantity) with dimension $(dimension(angle_quantity))")
    end
    return ustrip(uconvert(u"rad", angle_quantity))  # This returns Float64, dimensionless
end

# --- Time Conversions ---
function time_si_to_harmonic_units(time_si_quantity::Unitful.AbstractQuantity, harmonic_frequency)
    """
    Convert SI time to harmonic units (dimensionless).
    harmonic_frequency should be a Unitful quantity (e.g., MHz)
    Returns dimensionless time in units of 2π/ω
    """
    # Ensure input is a time quantity
    if dimension(time_si_quantity) != dimension(u"s")
        error("Input must be a time quantity. Got: $(time_si_quantity)")
    end
    
    # Convert harmonic frequency to angular frequency (rad/s)
    omega_si = 2π * ustrip(uconvert(u"Hz", harmonic_frequency))  # rad/s
    
    # Convert time to seconds then to harmonic units
    time_seconds = ustrip(uconvert(u"s", time_si_quantity))
    
    # t_harmonic = t_si * ω / (2π) = t_si * f (where f is frequency in Hz)
    return time_seconds * ustrip(uconvert(u"Hz", harmonic_frequency))
end

function time_harmonic_to_si(time_harmonic::Float64, harmonic_frequency, target_unit=u"s")
    """
    Convert harmonic time units to SI time.
    """
    # Convert to seconds first
    freq_hz = ustrip(uconvert(u"Hz", harmonic_frequency))
    time_seconds = time_harmonic / freq_hz
    
    # Convert to target unit
    return uconvert(target_unit, time_seconds * u"s")
end

function eigenvalues_au_to_harmonic_units(eigenvalues_au, harmonic_frequency)
    """
    Convert eigenvalues from atomic units to harmonic frequency units.
    This is the proper way to ensure energy-time consistency.
    
    IMPORTANT: The 2π factor is needed because our harmonic time is defined as
    t_harmonic = t_SI × f where f = ω/(2π) is frequency in Hz.
    """
    omega_au = frequency_to_energy_au_value(harmonic_frequency)
    # Include 2π factor for consistency with our time definition
    return eigenvalues_au ./ omega_au .* (2π)
end

end # module Units
