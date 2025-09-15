begin
	using Arrow
	using JLD2
	using Glob
	using Unitful
	using UnitfulAtomic
	using SparseArrays
	using Plots
	using LaTeXStrings
	using Latexify
	using UnitfulLatexify
	
	using FewMolecules
	using FewMolecules.Parameters
	using FewMolecules.Units
end

# Make sure you run this code from the root directory of the project

# # Define the electric dipole dimension if not already defined
# @derived_dimension ElectricDipole Unitful.𝐋*Unitful.𝐈*Unitful.𝐓

# Define the Debye unit (D) in terms of elementary charge and Bohr radius
# 1 D ≈ 0.393456 e⋅a₀ in atomic units as per the C++ implementation
# Note: e_au and a0_au are the elementary charge and Bohr radius in atomic units
@unit D "D" Debye 0.393456*UnitfulAtomic.e_au*UnitfulAtomic.a0_au false

# Register the unit to make it available to the u"" string macro
Unitful.register(Parameters)

begin
	load_params = FewMolecules.Parameters.Params()

	load_params.mode = 0 # Spectrum calculation mode
	load_params.vsWhat = 0 # 0: vs Electric field, 1: vs Magnetic field, 2: vs A factor, 3: vs Electric dipole, 4: vs Magnetic dipole
	# load_params.geom_values_flat = [] # One molecule: one at 0 nm, 0 degrees
	# load_params.geom_values_flat = [100.0, 60.0] # Two molecules: one at 0 nm, 0 degrees, one at 500 nm, 90 degrees
	# load_params.geom_values_flat = [100.0, 0.0] # Two molecules: one at 0 nm, 0 degrees, one at 500 nm, 90 degrees
	load_params.geom_values_flat = [500.0, 30.0, 500.0, 90.0] # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, 30 degrees, one at 500 nm, 90 degrees
	# load_params.geom_values_flat = [500.0, 90.0, 1000.0, 90.0] # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, 30 degrees, one at 500 nm, 90 degrees
	# load_params.geom_values_flat = [500.0, 0.0, 500.0, 180.0] # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, 30 degrees, one at 500 nm, 90 degrees
	# load_params.geom_values_flat = [100.0, 0.0, 100.0, 180.0] # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, 30 degrees, one at 500 nm, 90 degrees
	load_params.how_many_molecules = 3
	# load_params.E = 0.0u"kV/cm" # Electric field in kV/m
	# load_params.E = 11.48u"kV/cm" # Electric field in kV/m
	# load_params.E = 11.49u"kV/cm" # Electric field in kV/m
	# load_params.E = 11.45u"kV/cm" # Electric field in kV/m
	load_params.E = 11.44u"kV/cm" # Electric field in kV/m
	# load_params.E = 11.38u"kV/cm" # Electric field in kV/m
	# load_params.E = 11.40u"kV/cm" # Electric field in kV/m
	load_params.Phi_el = 0.0u"°"
	load_params.Theta_el = 90.0u"°"
	# load_params.B = 10.0u"Gauss"
	load_params.B = 50.0u"Gauss" # 5x the previous value, changed Zeeman definition
	# load_params.B = 0.0u"Gauss"
	load_params.Phi_mg = 0.0u"°"
	load_params.Theta_mg = 90.0u"°"
	load_params.A = 100.0u"kHz" # Spin-rotation coupling
	load_params.d_el = 0.1u"D" # Electric dipole
	load_params.d_mg = 10.0u"μB" # Magnetic dipole
	load_params.Mtot = -0.5 # Total angular momentum projection
	load_params.START_IT = 1
	# load_params.END_IT = 10001 # Single iteration
	load_params.END_IT = 501 # Single iteration
	# load_params.step = 0.00023 # Small step for quick test
	load_params.step = 0.00024 # Small step for quick test
	# load_params.step = 0.00046 # Small step for quick test
	# load_params.step = 0.015 # Small step for quick test
	# load_params.step = 0.0015 # Small step for quick test
	# load_params.step = 0.03 # Small step for quick test
	# load_params.step = 0.0002 # Small step for quick test
	# load_params.step = 0.00003 # Small step for quick test
	# load_params.step = 0.03 # Small step for quick test
	# load_params.step = 0.1 # Small step for quick test
	load_params.JMAX = 2 # Default from C++ constants.h
	load_params.spin = 5 # Default from C++ constants.h
	load_params.Brot = 0.29u"GHz" # Not used in this test
	load_params.hamiltonian_dtype = "Float64"
	load_params.eigenproblem_looped = true
	load_params.run_id_string = "S$(Int(load_params.spin))"
	load_params.print_colored = true;
	load_params.save_root_folder = ""
	load_params.initial_state_label = "Eigenstate 1365"
end 

gr()
# plot()

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, 
		framestyle=:box, 
		label=nothing, 
		grid=false,
		guidefont=16)
scalefontsizes(1.3)

field_to_latex_sig = Dict(
	"E" => L"\mathcal{E}",
	"B" => L"\mathcal{B}",
	"A" => L"\gamma",
	"d_el" => L"d_{el}",
	"d_mg" => L"d_{mg}",
)

# --- Conversion function using your Units module ---
function convert_energies_unitful(energies, from_unit::String, to_unit::String; Brot_value=0.29u"GHz")
	# Convert energies to atomic units (Hartree) first
	if from_unit == "atomic"
		energies_au = energies
	elseif from_unit == "Brot"
		Brot_au = Units.convert_param_to_au(Brot_value, :Brot)
		energies_au = energies .* ustrip(Brot_au)
	elseif from_unit == "Hz"
		hartree_to_hz = ustrip(auconvert(u"Hz", 1))
		energies_au = energies ./ hartree_to_hz
	elseif from_unit == "GHz"
		hartree_to_Ghz = ustrip(auconvert(u"GHz", 1))
		energies_au = energies ./ hartree_to_Ghz
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
	elseif to_unit == "GHz"
		hartree_to_Ghz = ustrip(auconvert(u"GHz", 1))
		energies_out = energies_au .* hartree_to_Ghz .*u"GHz"
	else
		error("Unsupported to_unit: $to_unit")
	end


	return energies_out
end;

function validate_loaded_file(fname, load_params)
	conds = [
		occursin(join(Int.(load_params.geom_values_flat), "_"), fname),
		occursin("Jmax$(load_params.JMAX)",fname),
		occursin("phiE$(round(Int, ustrip(load_params.Phi_el)))", fname),
		occursin("thetaE$(round(Int, ustrip(load_params.Theta_el)))", fname),
		occursin("phiB$(round(Int, ustrip(load_params.Phi_mg)))", fname),
		occursin("thetaB$(round(Int, ustrip(load_params.Theta_mg)))", fname),
		occursin("Mtot$(isinteger(load_params.Mtot) ? Int(round(load_params.Mtot, digits=0)) : load_params.Mtot)", fname),
	]
	if sum(conds) == length(conds)
		return true
	else
		# println(conds)
		return false
	end
end

function rectangle_from_coords(xb, xt, yb, yt)
	Shape(
		[xb, xt, xt, xb],
		[yb, yb, yt, yt],
	)
end

begin
	# Create subplot layout
	plot_spectrum = plot()
		
	data_folder = joinpath("./", FewMolecules.Utility.get_spectrum_filename(load_params, just_dir=true, override_params=true))
	# data_folder = joinpath("../", FewMolecules.Utility.get_spectrum_filename(load_params, just_dir=true, override_params=true))

	energies_folder = joinpath(data_folder, "energies")

	eigenstates_folder = joinpath(data_folder, "eigenstates")

	begin
		# List all energies file names to a List
		list_of_energy_files = Glob.glob(joinpath(energies_folder, "*.arrow"))
		list_of_eigenstate_files = Glob.glob(joinpath(eigenstates_folder, "*.jld2"))
		
		# Load all the values of the variable that was looped over
		global vsVal = FewMolecules.Utility.get_vs_folder_name(load_params.vsWhat)
	end

	if vsVal == "vsE"
		pattern = r"_E(\d+\.?\d*)_"
		unit_vs = u"kV/cm"
	elseif vsVal == "vsB"
		pattern = r"_B(\d+\.?\d*)_"
		unit_vs = u"Gauss"
	else
		@warn("This is not implemented yet")
	end

	geom_string = join(Int.(load_params.geom_values_flat), "_")


	begin
		# Extract E field values from filenames using regex "_E(\d+\.?\d*)_"
		valid_count = 0
		global energy_tables = []
		global vs_values = []
		for fname in list_of_energy_files
			m = match(pattern, fname)
			if m !== nothing
				if validate_loaded_file(fname, load_params)
					val = parse(Float64, m.captures[1]) 
					push!(vs_values, val * unit_vs)  # unit is u"kV/cm"
					push!(energy_tables, Arrow.Table(fname))
					valid_count += 1
				end
			else
				@warn("No $(vsVal[3:end]) field value found in filename: $fname")
			end
		end
		println("Number of valid files loaded: $valid_count")
		vs_values = uconvert.(unit_vs, vs_values)
		sorted_indices = sortperm(vs_values)
		global vs_values_sorted = vs_values[sorted_indices]
		global energy_tables_sorted = energy_tables[sorted_indices]
	end

	combined_energies = hcat([tbl.Eigenvalues for tbl in energy_tables_sorted]...)

	global from_unit = "Brot"  # Change if your data is in another unit

	# --- User input for Y axis unit ---
	global target_unit = "GHz"  # Options: "atomic", "Brot", "GHz, "Hz"

	unit_energies = target_unit == "atomic" ? u"hartree" :
					target_unit == "Brot"   ? unit(load_params.Brot) :
					target_unit == "GHz"   ? u"GHz" :
					target_unit == "Hz"     ? u"Hz" : nothing

	# --- Prepare data ---
	energies_converted = convert_energies_unitful(Matrix(combined_energies'), from_unit, target_unit)  # Transpose for plotting

	# Plot on left subplot (lines 277-279 style)
	for i in 1:size(energies_converted, 2)
		# lw_val = 1.2
		lw_val = 2
		# lw_val = 3
		alpha_val = 0.7
		plot!(plot_spectrum, vs_values_sorted, energies_converted[:, i], lw=lw_val, label=false, alpha=alpha_val)
	end
end


begin 
	labelfontsize=18
	tickfontsize=15
	annotsize=13
	# size_plot = (1200, 1550)
	# size_plot = (1200, 800)
	# size_plot = (1400, 1550)
	
	# labelfontsize=24
	# tickfontsize=20
	# annotsize=18
	# size_plot = (2400, 1600)
	# size_plot = (1600, 2000)
    # size_plot = (3000, 3000)
	size_plot = (1400, 1200)
	# size_plot = (2400, 3100)

	begin
		global xlabel = field_to_latex_sig[vsVal[3:end]] * " field (" * latexify(unit(vs_values[1])) * ")"
		global ylabel = target_unit == "atomic" ? "Energy (Hartree)" :
				target_unit == "Brot"   ? "Energy (Brot)" :
				target_unit == "GHz"   ? "Energy (GHz)" :
				target_unit == "Hz"     ? "Energy (Hz)" : "Energy"
	end

	begin 
		# Anticrossings for S=5/2
		previously_analyzed = [
			[(11.445, 11.54), (5.669, 5.673), "Anticrossing"], # Antycrossing 3 cząsteczki
			[(11.44,11.53), (10.821, 10.830), "Anticrossing"], # Anticrossing <--- 3 branche anticrossing.
			[(11.45,11.55), (6.5275, 6.5350), "Anticrossing"], # Antycrossing 3 cząsteczki
			[(11.425,11.575), (8.410,8.415), "Anticrossing"], # Antycrossing 3 cząsteczki
			[(11.44,11.55), (9.374,9.385), "Anticrossing"], # Anticrossing <--- 3 branche anticrossing.
			[(10.8, 10.85), (5.635, 5.645), "Crossing"], # Crossing 3 cząsteczki
			[(6.15, 6.25), (6.06, 6.09), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(7.55, 7.65), (8.365, 8.38), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(8.6, 8.7), (9.16, 9.17), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(8.99, 9.08), (9.515, 9.525), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(9.46, 9.58), (9.9725, 9.9825), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(8.52, 8.65), (10.6050, 10.612), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(7.55, 7.65), (12.54, 12.548), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(9.2, 9.3), (12.015, 12.045), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(9.65, 9.75), (12.42, 12.47), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(9.65, 9.75), (12.78, 12.8), "TBC (Likely crossing)"], # Raczej crossing
			[(7.87, 7.98), (12.915, 12.925), "TBC (Likely crossing)"], # Raczej crossing
			[(8.95, 9.1), (13.15, 13.162), "TBC (Anticrossing?)"], # Potencjalny antycrossing, policzyć dokładniej
			[(9.46, 9.58), (14.235, 14.246), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(8.965, 9.075), (14.595, 14.61), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(10.8, 10.85), (3.79, 3.8), "TBC (Anticrossing?)"], # Potencjalny antycrossing, policzyć dokładniej
			[(10.25, 10.34), (5.69, 5.705), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(14.4, 14.48), (4.6, 4.62), "TBC (Anticrossing?)"], # Potencjalny antycrossing, policzyć dokładniej
			[(10.25, 10.34), (5.69, 5.705), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(14.25,14.35), (6.94,6.98), "TBC (Likely crossing)"], # Raczej crossing
			[(13.76,13.84), (10.99,11.01), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej (ale raczej crossing)
			[(10.25,10.35), (10, 10.01), "TBC (Likely crossing)"], # Raczej crossing, policzyć dokładniej
			[(10.32,10.39), (9.926, 9.934), "TBC (Anticrossing?)"], # Potencjalny antycrossing, policzyć dokładniej
			[(12.1,12.2), (10.87, 10.882), "TBC (Anticrossing?)"], # Potencjalny antycrossing, policzyć dokładniej
			[(11.8, 11.95), (11.165, 11.18), "TBC (Anticrossing?)"], # Potencjalny antycrossing, policzyć dokładniej
		] 
		# ((11.42, 11.58), (9.445, 9.465))	<- Antycrossing Mtot=-1.5, S=5, before fix
	end 

	begin
		# global xlim, ylim = ((11.25, 11.70), (3.725,3.75)) 
		global xlim, ylim = ((11.45, 11.54), (10.330, 10.337)) # <---- 3 mol anticrossing, Mtot=-0.5
		# global xlim, ylim = ((11.44, 11.56), (7.469, 7.480)) # <---- 2 mol antictossing, 100 nm best
		# global xlim, ylim = ((11.4825, 11.503), (7.4743, 7.4757))  # <--- Zoomed in on the 2 mol antictossing
		# ylim = ustrip.(convert_energies_unitful.(ylim, "GHz", target_unit))
		# xlim, ylim, _ = previously_analyzed[1]
		# ylim = ustrip.(convert_energies_unitful.(ylim, "Brot", target_unit))
		# xlim, ylim = ((11.42, 11.58), (9.445, 9.465))	
		# xlim, ylim = ((0, 100), (-2.5, 14))	
		# xlim, ylim = ((10, 15), (5, 15))	
		# xlim, ylim = ((0, 15), (-6.5, 36.5))	
		# global xlim, ylim = (:auto, :auto)

		# MF spectrum
		# xlim, ylim = ((0, 100), (5, 20))	
	end 

	# # subset_coords = [((11.445, 11.54), (10.330, 10.337), "Zoom-in region")] # Zoom-in rectangle coordinates
	# subset_coords = [((11.4, 11.6), (10.275, 10.4), "Zoom-in region")] # Zoom-in rectangle coordinates

	# for entry in subset_coords
	# 	x_rng, y_rng, label = entry
	# 	if x_rng != xlim && y_rng != ylim
	# 		y_rng = ustrip.(convert_energies_unitful.(y_rng, "GHz", target_unit))
	# 		offset = ustrip.(convert_energies_unitful.(0.00005, "GHz", target_unit))
	# 	rect_x = x_rng[1]
	# 	rect_y = y_rng[1]
	# 	width = x_rng[2] - x_rng[1]
	# 	height = y_rng[2] - y_rng[1]
	# 	plot!(plot_spectrum, rectangle_from_coords(x_rng..., y_rng...), color="rgba(255, 255, 255, 0.01)", linecolor=:red, linewidth=2, label=false)
	# 	# annotate!(plot_spectrum, rect_x + width/2, rect_y + height + offset, (label, annotsize, :black, :center))
	# 	end	
	# end

	anticrossing_to_plot = 11.49914
	vline!([anticrossing_to_plot], lw=2, lc=:red, ls=:dash, label=L"\mathcal{E}^q=%$(anticrossing_to_plot)\,\mathrm{kV/cm}", alpha=0.66)

	outside_margin = 12Plots.mm

    # Format left subplot
    plot!(plot_spectrum,
        xlabel=xlabel, 
        ylabel=ylabel, 
        # title="(a)",
        title="",
        left_margin=outside_margin,
        right_margin=3*outside_margin,
        top_margin=outside_margin,
        bottom_margin=outside_margin,
        minorgrid=false,
        minorticks=10,
        minorgridalpha=0.1,
        minorgridlinewidth=1.0,
        minorgridstyle=:dash,
        ticks=:native,
        xlim=xlim, 
        ylim=ylim,
        grid=:all,
        gridstyle=:dash,
        gridalpha=0.25,
        gridlinewidth=1.5,
        framestyle=:box,
        framelinewidth=6,
        tickfontsize=labelfontsize,
        legendfontsize=annotsize,
		titlefontsize=1.5*labelfontsize,
        legendposition=:bottomleft,
		legend_background_color="rgba(255, 255, 255, 0.6)",
        ticklinewidth=4,
		size=size_plot, 
    )
    
    # Combine subplots
    final_plot = plot(plot_spectrum, layout=(1, 1), size=size_plot, dpi=600)
end

final_plot

begin 
# 	# savefig("plots/anticrossing_S8_vs_E$(join(xlim, "-"))$(unit_vs)_En$(join(ylim, "-"))$(target_unit).pdf")
# 	# savefig("plots/anticrossing_S8_vs_E$(join(xlim, "-"))$(unit_vs)_En$(join(ylim, "-"))$(target_unit).png")

	# savefig("plotting_scripts/plots/PDF/spectrum_highlighted_S$(load_params.spin)_$(vsVal)_singlemol_ThetaE$(load_params.Theta_el)_PhiE$(load_params.Phi_el)_ThetaB$(load_params.Theta_mg)_PhiB$(load_params.Phi_mg).pdf")
	# savefig("plotting_scripts/plots/PNG/spectrum_highlighted_S$(load_params.spin)_$(vsVal)_singlemol_ThetaE$(load_params.Theta_el)_PhiE$(load_params.Phi_el)_ThetaB$(load_params.Theta_mg)_PhiB$(load_params.Phi_mg).png")
	
	if xlim == :auto
		xlim = (minimum(vs_values_sorted), maximum(vs_values_sorted))
	end
	if ylim == :auto
		ylim = (minimum(energies_converted), maximum(energies_converted))
	end
	savefig(plot_spectrum, "plotting_scripts/plots/PDF/spectrum_S$(load_params.spin)_$(vsVal)$(join(xlim, "-"))$(unit_vs)_En$(join(ylim, "-"))$(target_unit)_ThetaE$(load_params.Theta_el)_PhiE$(load_params.Phi_el)_ThetaB$(load_params.Theta_mg)_PhiB$(load_params.Phi_mg)_$(join(Int.(load_params.geom_values_flat), "_")).pdf")
    savefig(plot_spectrum, "plotting_scripts/plots/PNG/spectrum_S$(load_params.spin)_$(vsVal)$(join(xlim, "-"))$(unit_vs)_En$(join(ylim, "-"))$(target_unit)_ThetaE$(load_params.Theta_el)_PhiE$(load_params.Phi_el)_ThetaB$(load_params.Theta_mg)_PhiB$(load_params.Phi_mg)_$(join(Int.(load_params.geom_values_flat), "_")).png")
end
 