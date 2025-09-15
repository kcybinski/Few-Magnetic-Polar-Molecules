using Arrow
using JLD2
using Glob
using Unitful
using UnitfulAtomic
using SparseArrays
# using GLMakie
using Plots

using FewMolecules
using FewMolecules.Parameters
using FewMolecules.Units

# Make sure you run this code from the root directory of the project

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

load_params = FewMolecules.Parameters.Params()

load_params.mode = 0 # Spectrum calculation mode
load_params.vsWhat = 0 # 0: vs Electric field, 1: vs Magnetic field, 2: vs A factor, 3: vs Electric dipole, 4: vs Magnetic dipole
# load_params.geom_values_flat = [] # One molecule: one at 0 nm, 0 degrees
# load_params.geom_values_flat = [500.0, 90.0] # Two molecules: one at 0 nm, 0 degrees, one at 500 nm, 90 degrees
load_params.geom_values_flat = [500.0, 30.0, 500.0, 90.0] # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, 30 degrees, one at 500 nm, 90 degrees
load_params.how_many_molecules = 3
# load_params.E = 0.0u"kV/cm" # Electric field in kV/m
load_params.E = 6.9u"kV/cm" # Electric field in kV/m
load_params.Phi_el = 0.0u"°"
load_params.Theta_el = 0.0u"°"
load_params.B = 10.0u"Gauss"
load_params.Phi_mg = 0.0u"°"
load_params.Theta_mg = 0.0u"°"
load_params.A = 100.0u"kHz" # Spin-rotation coupling
load_params.d_el = 0.1u"D" # Electric dipole
load_params.d_mg = 10.0u"μB" # Magnetic dipole
load_params.Mtot = -1 # Total angular momentum projection
load_params.START_IT = 1
load_params.END_IT = 101 # Single iteration
load_params.step = 0.003 # Small step for quick test
# load_params.START_IT = 1
# load_params.END_IT = 11 # Single iteration
# load_params.step = 0.01 # Small step for quick test
load_params.JMAX = 2 # Default from C++ constants.h
load_params.spin = 5 # Default from C++ constants.h
load_params.Brot = 0.29u"GHz" # Not used in this test
load_params.hamiltonian_dtype = "Float64"
load_params.eigenproblem_looped = true
load_params.run_id_string = "s$(Int(load_params.spin))"
load_params.print_colored = true;

data_folder = joinpath("./Julia_Few_Molecules", FewMolecules.Utility.get_spectrum_filename(load_params, just_dir=true, override_params=true))

energies_folder = joinpath(data_folder, "energies")
eigenstates_folder = joinpath(data_folder, "eigenstates")

simulation_params = load(joinpath(data_folder, "simulation_params.jld2"))["params"]


# List all energies file names to a List
list_of_energy_files = Glob.glob(joinpath(energies_folder, "*.arrow"))
list_of_eigenstate_files = Glob.glob(joinpath(eigenstates_folder, "*.jld2"))

# Load all the values of the variable that was looped over
vsVal = FewMolecules.Utility.get_vs_folder_name(load_params.mode)

if vsVal == "vsE"
    pattern = r"_E(\d+\.?\d*)_"
    unit_vs = u"kV/cm"
elseif vsVal == "vsB"
    pattern = r"_B(\d+\.?\d*)_"
    unit_vs = u"Gauss"
else
    @warn("This is not implemented yet")
end

# Extract E field values from filenames using regex "_E(\d+\.?\d*)_"
energy_tables = []
vs_values = []
for fname in list_of_energy_files
    m = match(pattern, fname)
    if m !== nothing
        val = parse(Float64, m.captures[1]) 
        push!(vs_values, val * unit_vs)  # unit is u"kV/cm"
        push!(energy_tables, Arrow.Table(fname))
    else
        @warn("No $(vsVal[3:end]) field value found in filename: $fname")
    end
end
vs_values = uconvert.(unit_vs, vs_values)
sorted_indices = sortperm(vs_values)
vs_values_sorted = vs_values[sorted_indices]
energy_tables_sorted = energy_tables[sorted_indices]

combined_energies = hcat([tbl.Eigenvalues for tbl in energy_tables_sorted]...);

# --- Conversion function using your Units module ---
function convert_energies_unitful(energies::Matrix{<:Real}, from_unit::String, to_unit::String; Brot_value=0.29u"GHz")
    # Convert energies to atomic units (Hartree) first
    if from_unit == "atomic"
        energies_au = energies
    elseif from_unit == "Brot"
        Brot_au = Units.convert_param_to_au(Brot_value, :Brot)
        energies_au = energies .* ustrip(Brot_au)
    elseif from_unit == "Hz"
        hartree_to_hz = ustrip(auconvert(u"Hz", 1))
        energies_au = energies ./ hartree_to_hz
    else
        error("Unsupported from_unit: $from_unit")
    end

    # Now convert atomic units to target unit
    if to_unit == "atomic"
        energies_out = energies_au.*u"hartree"
    elseif to_unit == "Brot"
        Brot_au = Units.convert_param_to_au(Brot_value, :Brot)
        energies_out = energies_au ./ ustrip(Brot_au) .* unit(Brot_value)
    elseif to_unit == "Hz"
        hartree_to_hz = ustrip(auconvert(u"Hz", 1))
        energies_out = energies_au .* hartree_to_hz .*u"Hz"
    else
        error("Unsupported to_unit: $to_unit")
    end


    return energies_out
end;

# --- User input for Y axis unit ---
target_unit = "Brot"  # Options: "atomic", "Brot", "Hz"
unit_energies = target_unit == "atomic" ? u"hartree" :
                target_unit == "Brot"   ? unit(load_params.Brot) :
                target_unit == "Hz"     ? u"Hz" : nothing

# --- Prepare data ---
from_unit = "Brot"  # Change if your data is in another unit
energies_converted = convert_energies_unitful(Matrix(combined_energies'), from_unit, target_unit)  # Transpose for plotting


# --- Plot ---
gr()
plot()
for i in 1:size(energies_converted, 2)
    plot!(vs_values_sorted, energies_converted[:, i], lw=2, label=false)
    # plot!(vs_values_sorted, energies_converted[:, i], lw=0, marker=:circle, markersize=2, label=false)
end

xlabel = string("$(vsVal[3:end]) Field(", unit(vs_values[1]), ")")
ylabel = target_unit == "atomic" ? "Energy (Hartree)" :
         target_unit == "Brot"   ? "Energy (Brot)" :
         target_unit == "Hz"     ? "Energy (Hz)" : "Energy"

# xlim_plt1 = (6.875, 7.175)
xlim_plt1 = (7.03, 7.08)
# xlim_plt1 = (11.425 ,11.575)
# ylim_plt1 = (5.669,5.673)
ylim_plt1 = (9.607, 9.616)

# size_plot = (3000, 2400)
size_plot = (1200,800)
# plot!(xlabel=xlabel, ylabel=ylabel, legend=false, xlim=xlim, ylim=ylim, size=size_plot)
plot!(
    xlabel=xlabel, 
    ylabel=ylabel, 
    legend=false, 
    size=size_plot,
    minorgrid=true,
    minorticks=10,
    ticks=:native,
    xlim=xlim_plt1, 
    ylim=ylim_plt1,
    tickfontsize=15,
    grid=:all,
    gridstyle=:dash,
    gridlinewidth=2,
    )

display(current())

xlim = (11.45, 11.55)
ylim = (5.669,5.673)

# Load channels
channels = Arrow.Table(joinpath(data_folder, "channels.arrow"))

vs_vals_in_xlim = ((xlim[1]*unit_vs) .< vs_values_sorted .< (xlim[2]*unit_vs))
println("Entries in X lim: $(sum(vs_vals_in_xlim))")

energies_in_ylim = ((ylim[1]*unit_energies) .< energies_converted[vs_vals_in_xlim, :] .< (ylim[2]*unit_energies));
println("Energies in Y lim : $(sum(energies_in_ylim)), size: $(size(energies_in_ylim))")

eigenstate_inds = spzeros(Int, size(energies_in_ylim)[1], size(energies_in_ylim)[2])
for i in 1:size(eigenstate_inds)[1]
    for j in 1:size(eigenstate_inds)[2]
        if energies_in_ylim[i,j] > 0
            eigenstate_inds[i,j] = Int(j)
        end
    end
end
energies_window_min = minimum(eigenstate_inds[eigenstate_inds .> 0])
energies_window_max = maximum(eigenstate_inds[eigenstate_inds .> 0])

if vsVal == "vsE"
    pattern_eigenstates = r"_E(\d+\.?\d*)_"
elseif vsVal == "vsB"
    pattern_eigenstates = r"_B(\d+\.?\d*)_"
else
    @warn("This is not implemented yet")
end

es_array = []
for fname in list_of_eigenstate_files
    m = match(pattern, fname)
    if m !== nothing
        val = parse(Float64, m.captures[1]) 
        if val*unit_vs in vs_values_sorted[vs_vals_in_xlim]
            eigenstates = jldopen(fname)["vectors"][:, energies_window_min:energies_window_max]
            push!(es_array, eigenstates)
            # println("Found for val=$(val), shape: $(size(eigenstates))")
        end
    else
        continue
        # @warn "No $(vsVal[3:end]) field value found in filename: $fname"
    end
end

function build_eigenstate_hover_text(eigenstates, energies, vs_value_fun, basis, es_ind, es_ind_min;threshold_ratio=0.1, join_lines_str="<br>")
    amps = abs2.(eigenstates)
    es_energy = energies[es_ind_min+es_ind-1]
    max_amp = maximum(amps)
    threshold = 0.1 * max_amp

    amps[amps .< threshold] .= 0

    dropzeros!(amps)

    coords = (amps.nzind, amps.nzval)

    sortinds = sortperm(coords[2], rev=true)  # Sort by amplitude in descending order
    coords = (coords[1][sortinds], coords[2][sortinds])
    
    lines = String[]
    push!(lines, "Eigenstate no. $(es_ind_min+es_ind-1)")
    push!(lines, "vsVal: $(vs_value_fun), Energy: $(round(ustrip(es_energy), digits=3))")
    for (ind, amp) in zip(coords[1], coords[2])
        push!(lines, "|Ψ|²=$(round(amp, digits=3)), basis=$(basis[ind])")
    end
    
    return(join(lines, join_lines_str))
end

hover_matrix = []
for (pt_ind, es_arr) in enumerate(es_array)
    col_tab = []
    energies = energies_converted[vs_vals_in_xlim, :][pt_ind, :]
    vs_value = vs_values_sorted[vs_vals_in_xlim][pt_ind]
    for (es_ind, es_col) in enumerate(eachcol(es_arr))
        es_col = SparseVector(es_col)
        push!(col_tab, build_eigenstate_hover_text(es_col, energies, vs_value, channels, es_ind, energies_window_min))
    end
    push!(hover_matrix, Array(col_tab))
end
hover_matrix = reduce(hcat, hover_matrix)

# Plotting the subset with hover text
plotlyjs()
plot()
for i in 1:size(energies_converted[vs_vals_in_xlim, energies_window_min:energies_window_max], 2)
    plot!(
        vs_values_sorted[vs_vals_in_xlim], 
        energies_converted[vs_vals_in_xlim, energies_window_min:energies_window_max][:, i], 
        lw=2, 
        label=false,
        hover=hover_matrix[i]
        )
    # plot!(vs_values_sorted, energies_converted[:, i], lw=0, marker=:circle, markersize=2, label=false)
end

xlabel = string("$(vsVal[3:end]) Field(", unit(vs_values[1]), ")")
ylabel = target_unit == "atomic" ? "Energy (Hartree)" :
         target_unit == "Brot"   ? "Energy (Brot)" :
         target_unit == "Hz"     ? "Energy (Hz)" : "Energy"

# size_plot = (3000, 2400)
size_plot = (1000,800)
plot!(xlabel=xlabel, ylabel=ylabel, legend=false, xlim=xlim, ylim=ylim, size=size_plot)
# plot!(xlabel=xlabel, ylabel=ylabel, legend=false, size=size_plot)
display(current())

