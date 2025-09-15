module PhysicalConstantsJL

using PhysicalConstants
using Unitful # Add Unitful
using UnitfulAtomic # Add UnitfulAtomic

const Constants = PhysicalConstants.CODATA2022

# Export order doesn't strictly matter for functionality but can be organized for readability.
# These are already Unitful.Quantity objects from PhysicalConstants.CODATA2022
export fine_structure_constant, bohr_radius, elementary_charge, electron_mass
export speed_of_light_in_vacuum, planck_constant, boltzmann_constant
export vacuum_mag_permeability, vacuum_elec_permittivity, hartree_energy

# Fundamental Physical Constants from CODATA2022, in SI units (and are Unitful Quantities)

# Electron mass (kg)
const electron_mass = Constants.m_e

# Fine-structure constant (dimensionless)
# Ensure it's treated as a pure number if it's not already, or handle appropriately.
# PhysicalConstants.CODATA2022.α is already a Float64, which is fine for a dimensionless constant.
const fine_structure_constant = Constants.α

# Speed of light in vacuum (meters/second)
const speed_of_light_in_vacuum = Constants.c_0

# Hartree energy (Joules)
# Defined using UnitfulAtomic.hartree and converted to Joules for consistency with other SI unit constants in this file.
const hartree_energy = uconvert(u"J", 1u"hartree")

# Bohr radius (meters)
# PhysicalConstants.CODATA2018.a_0 is available.
# UnitfulAtomic.bohr is also available.
const bohr_radius = Constants.a_0

# Elementary charge (Coulombs)
const elementary_charge = Constants.e

# Planck constant (Joule-second)
const planck_constant = Constants.h # This is h
# For ħ (h-bar), use Constants.ħ or planck_constant / (2π)
export ħ
const ħ = Constants.ħ # Reduced Planck constant

# Boltzmann constant (Joules/Kelvin)
const boltzmann_constant = Constants.k_B

# Vacuum magnetic permeability (N/A^2 or H/m)
const vacuum_mag_permeability = Constants.μ_0

# Vacuum electric permittivity (F/m)
const vacuum_elec_permittivity = Constants.ε_0

# Atomic units (can be explicitly imported or used via UnitfulAtomic)
# Example: export au_energy, au_length
# These are typically Unitful.hartree and Unitful.bohr respectively if UnitfulAtomic is used,
# or defined via UnitfulAtomic.au"Hartree" etc.

# Other useful constants can be added here as needed, for example:
# Avogadro constant (1/mol)
# const avogadro_constant = N_A

# Molar gas constant (J/(mol*K))
# const molar_gas_constant = R

end # module PhysicalConstantsJL
