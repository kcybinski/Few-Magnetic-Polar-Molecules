begin
	using Arrow
	using JLD2
	using Glob
	using Unitful
	using UnitfulAtomic
	using SparseArrays
	using Plots
	using PlutoUI
	using CSV
	using DataFrames
	using Interpolations
	using FFTW
	using LaTeXStrings
	using Latexify
	using UnitfulLatexify

	using FewMolecules
	using FewMolecules.Parameters
	using FewMolecules.Units
	using FewMolecules.Utility
end

# Make sure you run this code from the root directory of the project

function plot_fft_for_observable(obs_name, ω, dfs_observables, plot_dir, save_plots=true; left_labels=false,maxfreq=nothing)

	if occursin("mtot", obs_name)
		mol = obs_name[end]
		m_label = "m$(mol)"
		ms_label = "ms$(mol)"
		t = dfs_observables[m_label].time_harmonic
		y = dfs_observables[m_label].value + dfs_observables[ms_label].value
	else
		t = dfs_observables[obs_name].time_harmonic
		y = dfs_observables[obs_name].value
	end

    # Compute the DFT
    Y = fft(y)

    # Compute frequency axis properly scaled to trap frequency
    N = length(t)
    tmax_harmonic = t[end]     # max time in harmonic units

    # Print time step information for clarity
    dt_harmonic = t[2] - t[1]  # time step in harmonic units
    max_freq = 1 / (2 * dt_harmonic)  # Nyquist frequency in harmonic units
    println("Time step: $(dt_harmonic) × (2π/ω)")
    println("This limits detectable frequencies to $(1/(2*dt_harmonic))ω")
    println("Trap frequency: 2π × $(ustrip(ω)/(2π)) kHz, period: $(ustrip(2π/ω)*1000) ms/(2π/ω)")

    # In harmonic units, frequency will be directly in units of ω
    # freq_harmonic = (0:N-1) ./ (tmax_harmonic + dt_harmonic) 
    # For a time series in harmonic units, frequency should be expressed in units of ω
    freq_harmonic = (0:N-1) ./ (N * dt_harmonic)  # Frequencies in units of ω

    # Only show the positive frequencies up to Nyquist frequency
    max_idx = div(N, 2) + 1
    plot_freq = freq_harmonic[1:max_idx]
    plot_amp = 2 .* abs.(Y[1:max_idx]) ./ N  # Normalize by N and multiply by 2 (except DC)
    plot_amp[1] = plot_amp[1] / 2  # Don't double-count DC component

	min_amp = floor(log10(minimum(plot_amp[plot_amp .> 0])))  # Use floor instead of round
	max_amp = ceil(log10(maximum(plot_amp)))                   # Use ceil instead of round

	# Generate full range of integer powers
	full_range = Int(min_amp):Int(max_amp)                     # Full possible range

	if length(full_range) <= 3
		xticks_powers = full_range                             # Use all if 3 or fewer
	else
		# Choose min, max, and closest to midpoint for maximum coverage
		min_power = full_range[1]                              # Minimum power
		max_power = full_range[end]                            # Maximum power
		mid_target = (min_power + max_power) / 2               # Ideal midpoint
		mid_power = full_range[argmin(abs.(full_range .- mid_target))]  # Closest to midpoint
		
		xticks_powers = [min_power, mid_power, max_power]      # Three powers with max coverage
	end

	# Create all powers for minor grid
	all_powers = Int(min_amp):Int(max_amp)
	all_ticks = [10.0^p for p in all_powers]

	# Create main ticks and labels
	xticks = [10.0^p for p in xticks_powers]                  # Create ticks as powers of 10
	xticks_labels = [L"10^{%$(p)}" for p in xticks_powers]    # Create LaTeX labels


	p = plot(
		ylabel="",
		xlabel="Amplitude", 
		xscale=:log10,
		framestyle=:box,
        xmirror=true,
		xticks=(xticks, xticks_labels), 
		# Add minor ticks at all powers for grid
		minorticks=all_ticks,
		minorgrid=false,
		minorgridalpha=0.3,
		xrot=0,
		ymirror=!left_labels ? true : false,
		bottom_margin=save_plots ? 2Plots.mm : 1Plots.mm,
        left_margin=save_plots ? 2Plots.mm : left_labels ? 15Plots.mm : 1Plots.mm,
        top_margin=save_plots ? 15Plots.mm : 5Plots.mm,
        right_margin=save_plots ? 15Plots.mm : left_labels ? 1Plots.mm : 15Plots.mm,
		)
		
	# Plot as main plot
	plot!(
		p,
		plot_amp,
		plot_freq,
		ylabel="Frequency " * L"(\omega_{\mathrm{trap}})", 
		legend=false,
		xscale=:log10,
		# xticks=xticks,
		xlim=(10.0^(min_amp - 0.5), 10.0^(max_amp + 0.5)),
		framestyle=:box,
		title=save_plots ? "FFT of $obs_name" : "",
		ylim=(-0.01, maxfreq !== nothing ? maxfreq : max_freq),
	)

	if !save_plots
		return p
	end
    savefig(joinpath(joinpath(plot_dir, "png"), "quench_fft_$(obs_name).png"))
    savefig(joinpath(joinpath(plot_dir, "pdf"), "quench_fft_$(obs_name).pdf"))
end

# Define the Debye unit (D) in terms of elementary charge and Bohr radius
# 1 D ≈ 0.393456 e⋅a₀ in atomic units as per the C++ implementation
# Note: e_au and a0_au are the elementary charge and Bohr radius in atomic units
@unit D "D" Debye 0.393456*UnitfulAtomic.e_au*UnitfulAtomic.a0_au false

# Register the unit to make it available to the u"" string macro
Unitful.register(Parameters)

begin
	load_params = FewMolecules.Parameters.Params()

	load_params.mode = 1 # Quench calculation mode
	load_params.vsWhat = 1 # 0: vs Electric field, 1: vs Magnetic field, 2: vs A factor, 3: vs Electric dipole, 4: vs Magnetic dipole
	load_params.geom_values_flat = [] # One molecule: one at 0 nm, 0 degrees
	# load_params.geom_values_flat = [100.0, 60.0] # Two molecules: one at 0 nm, 0 degrees, one at 500 nm, 90 degrees
	# load_params.geom_values_flat = [100.0, 0.0, 100.0, 180.0] # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, 30 degrees, one at 500 nm, 90 degrees
	# load_params.geom_values_flat = [500.0, 90.0, 1000.0, 90.0] # Three molecules: one at 0 nm, 0 degrees, one at 500 nm, 90 degrees, one at 1000 nm, 90 degrees
	load_params.how_many_molecules = length(load_params.geom_values_flat) ÷ 2 + 1
	load_params.E = 11.4848u"kV/cm" # Electric field in kV/m
	load_params.Phi_el = 0.0u"°"
	load_params.Theta_el = 90.0u"°"
	load_params.B = 50.0u"Gauss"
	load_params.Phi_mg = 0.0u"°"
	load_params.Theta_mg = 90.0u"°"
	load_params.A = 100.0u"kHz" # Spin-rotation coupling
	load_params.d_el = 0.1u"D" # Electric dipole
	load_params.d_mg = 10.0u"μB" # Magnetic dipole
	
	load_params.JMAX = 2 # Default from C++ constants.h
	load_params.spin = 5 # Default from C++ constants.h
	load_params.Brot = 0.29u"GHz" # Not used in this test
	load_params.hamiltonian_dtype = "ComplexF64"

	load_params.harmonic_frequency = 100 * 2π # in kHz
	load_params.time_is_harmonic = true
	load_params.tmax = 500.0
	load_params.dt = 0.25
	# load_params.ITS_PER_OVERRIDE = 6


	plot_in_ms = true
	time_in_ms = load_params.time_is_harmonic ? load_params.tmax * u"ms" * (2π / load_params.harmonic_frequency) : load_params.tmax * u"ms"
	time_tab_in_ms = collect(0.0:load_params.dt:load_params.tmax) .* (load_params.time_is_harmonic ? u"ms" * (2π / load_params.harmonic_frequency) : u"ms")
	
	# load_params.quenched_E = 11.4345u"kV/cm"
	# load_params.quenched_E = 11.4927u"kV/cm"
	load_params.quenched_E = 11.49914u"kV/cm"
	load_params.quenched_R = [(100.0u"nm", 0u"°")]
	# load_params.quenched_R = [(100.0u"nm", 0.0u"°"), (100.0u"nm", 180.0u"°")]

	angle_var = load_params.quenched_R[1][2] # Extract angle for labeling
	lab_var = L"(a)"
	# load_params.observables = ["M", "Ms", "m1", "ms1", "m2", "ms2", "m3", "ms3"]
	# load_params.observables = ["M", "Ms", "m1", "ms1", "m2", "ms2"]
	load_params.observables = ["M", "Ms", "m1", "ms1"]
	# load_params.initial_eigenstate = 1884
	# load_params.initial_state_label = "Eigenstate $(load_params.initial_eigenstate)"


	# 1 molecule breakdown of the interesting quench
	# load_params.initial_state_label = "Product State [1.0;-1.0;2.5;-2.5]" # [j¹ₘₐₓ, mⱼ¹, s¹, mₛ¹]
	# load_params.Mtot = -3.5 # Total angular momentum projection

	# load_params.initial_state_label = "Product State [1.0;1.0;2.5;2.5]" # [j¹ₘₐₓ, mⱼ¹, s¹, mₛ¹]
	# load_params.Mtot = 3.5 # Total angular momentum projection

	load_params.initial_state_label = "Product State [1.0;-1.0;2.5;0.5]" # [j¹ₘₐₓ, mⱼ¹, s¹, mₛ¹]
	load_params.Mtot = -0.5 # Total angular momentum projection


	# 2 molecule breakdowns of the interesting quench
	# load_params.initial_state_label = "Product State [1.0;-1.0;2.5;-2.5;1.0;1.0;2.5;2.5]"
	# load_params.Mtot = 0.0 # Total angular momentum projection
	
	# load_params.initial_state_label = "Product State [1.0;-1.0;2.5;-2.5;1.0;-1.0;2.5;0.5]"
	# load_params.Mtot = -4.0 # Total angular momentum projection
	
	# load_params.initial_state_label = "Product State [1.0;1.0;2.5;2.5;1.0;-1.0;2.5;0.5]"
	# load_params.Mtot = 3.0 # Total angular momentum projection


	# 3 molecule interesting quench
	load_params.initial_state_label = "Product State [1.0;-1.0;2.5;-2.5;1.0;1.0;2.5;2.5;1.0;-1.0;2.5;0.5]"
	load_params.Mtot = -0.5 # Total angular momentum projection


	# load_params.productstate = []
	load_params.run_id_string = "S$(Int(load_params.spin))"
	load_params.print_colored = true;
	load_params.save_root_folder = ""

	counter_double_counting = false

	if load_params.ITS_PER_OVERRIDE !== nothing
		num_iterations = Int(rationalize(load_params.tmax / (load_params.dt * load_params.ITS_PER_OVERRIDE)))
	else
		num_iterations = Int(rationalize(load_params.tmax / load_params.dt)) + 1
	end
end

begin
	gr()

	plot_font = "Computer Modern"
	default(fontfamily=plot_font,
			linewidth=1, 
			framestyle=:box, 
			label=nothing, 
			grid=false,
			guidefont=14,
			tickfont=14,
			legendfont=14,
			titlefont=16,
			size=(1200, 800),
			left_margin=15Plots.mm,
			right_margin=15Plots.mm,
			top_margin=8Plots.mm,
			bottom_margin=8Plots.mm
			)
	scalefontsizes(1.3)
	# xlim = (0, 25) # <--- clear plotting of dense dynamics
	# xlim = (0, 100) # <--- clear plotting of dense dynamics
	# xlim = (0, 120) # <--- 3 molecules old dynamics
	xlim = (0, 250) # <--- 3 molecules new dynamics # Adjusted x-axis limits for consistency
	# xlim = (0, 5000) # <-- 2 molecules long dynamics
	# xlim = (0, 500) # <-- 2 molecules regular dynamics

	maxfreq_fft = 1.0
end

obs_rename_dict = Dict(
	"M" => L"M_J",
	"Ms" => L"M_S",
	"m1" => L"m_{j}^{(1)}",
	"ms1" => L"m_{s}^{(1)}",
	"m2" => L"m_{j}^{(2)}",
	"ms2" => L"m_{s}^{(2)}",
	"m3" => L"m_{j}^{(3)}",
	"ms3" => L"m_{s}^{(3)}",
)

begin
	data_folder = joinpath("./",FewMolecules.Utility.get_quench_save_folder(load_params, override=true))
	# data_folder = "./data/quench_results/1_molecule/EF11/MF10/A100/D_EL0.1/D_MG10/Eq11.493_Rq/Eigenstate $(load_params.initial_eigenstate)/tmax_100.0dt_0.1unit_harmonic"

	# # List all observable files
	# list_of_observable_files = Glob.glob(joinpath(data_folder, "*.csv"))
	list_of_observable_files = [f for f in readdir(data_folder, join=true) if occursin(".csv", f)]

	dfs_observables = Dict([
		(obs => DataFrame()) for obs in load_params.observables
	])

	channels = Arrow.Table(joinpath(data_folder, "channels_S$(load_params.spin)_Jmax$(load_params.JMAX)_Mtot$(load_params.Mtot).arrow"))

	for filename_tmp in list_of_observable_files
		if occursin("couplings", filename_tmp)
			println("   Skipping coupling file: $filename_tmp")
			continue  # Skip coupling files
		end

		for obs in load_params.observables
			try
				dfs_observables[obs] = CSV.read(joinpath(data_folder, "quench_observable_$obs.csv"), DataFrame; comment="#")
			catch e
				try 
					dfs_observables[obs] = CSV.read(joinpath(data_folder, "quench_$obs.csv"), DataFrame; comment="#")
				catch e_inner
					println("   Error reading $obs in second attempt: $e_inner")
				end
				println("   Error reading $obs: $e")
			end
		end
	end

	plot_dir = joinpath(data_folder, "plots")
	isdir(plot_dir) ? nothing : mkdir(plot_dir)
	plot_png_dir = joinpath(plot_dir, "PNG")
	isdir(plot_png_dir) ? nothing : mkdir(plot_png_dir)
	plot_pdf_dir = joinpath(plot_dir, "PDF")
	isdir(plot_pdf_dir) ? nothing : mkdir(plot_pdf_dir)

	gr()

	for obs in load_params.observables
		plot(
			dfs_observables[obs].time_harmonic,
			dfs_observables[obs].value,
			title="Observables vs Time",
			label=obs,
		)
		plot!(
		title="Observables vs Time",
		xlabel="Time " * L"(2\pi/\omega_{\mathrm{trap}})",
		ylabel=L"\langle %$(obs_rename_dict[obs]) \rangle",
		legend=true,
		xlim=xlim,
		)
		savefig(joinpath(plot_png_dir, "quench_observable_$(obs).png"))
		savefig(joinpath(plot_pdf_dir, "quench_observable_$(obs).pdf"))

		plot_fft_for_observable(obs, load_params.harmonic_frequency, dfs_observables, plot_dir, true)
	end

	begin
		plot()
		for obs in load_params.observables
			plot!(
				dfs_observables[obs].time_harmonic,
				dfs_observables[obs].value,
				title="Observables vs Time",
				label=L"%$(obs_rename_dict[obs])",
				palette=:tab10,
			)
		end

		plot!(
			dfs_observables["M"].time_harmonic,
			dfs_observables["M"].value + dfs_observables["Ms"].value,
			label=L"M_J + M_S",
			linestyle=:dash
		)

		plot!(
			title="Observables vs Time",
			xlabel="Time " * L"(2\pi/\omega_{\mathrm{trap}})",
			ylabel=L"\langle \Psi(t) | \hat{O}~| \Psi(t) \rangle",
			legend=true,
			xlim=xlim,
			legend_background_color="rgba(255, 255, 255, 0.6)",
		)

		savefig(joinpath(plot_png_dir, "00_quench_observables_all.png"))
		savefig(joinpath(plot_pdf_dir, "00_quench_observables_all.pdf"))
	end
end 

begin
	# Plot the total angular momentum projection M_J and total spin projection M_S on the same plot with two y-axes
	# color_1 = :cyan3
	# color_2 = :magenta4
	color_1 = palette(:tab10)[1]
	color_2 = palette(:tab10)[2]
	lw = 1.3
	xlim_tmp = (0, 50)

	p_M_Ms = plot(
		dfs_observables["M"].time_harmonic,
		dfs_observables["M"].value,
		xlim=xlim_tmp,
		label=L"M_J",
		ylabel=L"\langle M_J \,\rangle",
		xlabel="Time " * L"(2\pi/\omega_{\mathrm{trap}})",
		lw=lw,
		color=color_1,
		y_guidefontcolor=color_1,
		y_foreground_color_axis=color_1,
		y_foreground_color_text=color_1,
		y_foreground_color_border=color_1,
	)
	

	plot!(
		twinx(),
		dfs_observables["Ms"].time_harmonic,
		dfs_observables["Ms"].value,
		xlim=xlim_tmp,
		label=L"M_S",
		ylabel=L"\langle M_S \rangle",
		lw=lw,
		color=color_2,
		y_guidefontcolor=color_2,
		y_foreground_color_axis=color_2,
		y_foreground_color_text=color_2,
		y_foreground_color_border=color_2,
	)

	
	plot!(
		legend=false,
		xlim=xlim_tmp,
		frame=:box,
	)
		
	plot_M_fft = plot_fft_for_observable("M", load_params.harmonic_frequency, dfs_observables, plot_dir, false;left_labels=true, maxfreq=maxfreq_fft)
	plot_Ms_fft = plot_fft_for_observable("Ms", load_params.harmonic_frequency, dfs_observables, plot_dir, false; maxfreq=maxfreq_fft)
	
	plot!(
		plot_M_fft,
		xmirror=false,
		bottom_margin=16Plots.mm, # <--- TODO Change to 16Plots.mm if single molecule plotting
		left_margin=15Plots.mm,
	)
	plot!(
		plot_Ms_fft,
		xmirror=false,
		bottom_margin=16Plots.mm,
		right_margin=15Plots.mm,
	)

	final_plot = plot(
		stack(([plot_M_fft], [p_M_Ms], [plot_Ms_fft]), dims=1)...,
		layout=grid(1, 3, widths=(1/6, 4/6, 1/6)),
		legend_background_color="rgba(255, 255, 255, 0.6)",
		size=(2000, 900),
	)

	savefig(joinpath(plot_png_dir, "04_quench_observables_M_Ms.png"))
	savefig(joinpath(plot_pdf_dir, "04_quench_observables_M_Ms.pdf"))
end

begin
	plot_vec_tmp = []
	fft_mj_vec_tmp = []
	fft_ms_vec_tmp = []
	for mol in 1:load_params.how_many_molecules
		# color_1 = :cyan3
		# color_2 = :magenta4
		color_1 = palette(:tab10)[1]
		color_2 = palette(:tab10)[2]
		lw = 1.3

		m_label = "m$(mol)"
		ms_label = "ms$(mol)"

		p_tmp = plot(
			dfs_observables[m_label].time_harmonic,
			dfs_observables[m_label].value,
			label=L"m_{j}^{(%$(mol))}",
			ylabel=L"\langle m_{j}^{(%$(mol))} \,\rangle",
			xlabel="Time " * L"(2\pi/\omega_{\mathrm{trap}})",
			lw=lw,
			color=color_1,
			y_guidefontcolor=color_1,
			y_foreground_color_axis=color_1,
			y_foreground_color_text=color_1,
			y_foreground_color_border=color_1,
		)

		plot!(
			twinx(),
			dfs_observables[ms_label].time_harmonic,
			# dfs_observables[ms_label].value,
			ms_label == "ms1" ? dfs_observables[ms_label].value .- 9.999999e-1 : dfs_observables[ms_label].value,
			# ylabel=L"\langle m_{s}^{(%$(mol))} \rangle",
			ylabel=ms_label == "ms1" ? L"\langle m_{s}^{(%$(mol))}\rangle - 9.999999\cdot10^{-1}" : L"\langle m_{s}^{(%$(mol))} \rangle",
			lw=lw,
			color=color_2,
			y_guidefontcolor=color_2,
			y_foreground_color_axis=color_2,
			y_foreground_color_text=color_2,
			y_foreground_color_border=color_2,
		)

		plot!(
			legend=false,
			xlim=(0, 50),
			frame=:box,
		)

		savefig(joinpath(plot_png_dir, "quench_observables_m$(mol)_ms$(mol).png"))
		savefig(joinpath(plot_pdf_dir, "quench_observables_m$(mol)_ms$(mol).pdf"))

		plot_mj_fft_tmp = plot_fft_for_observable(m_label, load_params.harmonic_frequency, dfs_observables, plot_dir, false;left_labels=true, maxfreq=maxfreq_fft)
		plot_ms_fft_tmp = plot_fft_for_observable(ms_label, load_params.harmonic_frequency, dfs_observables, plot_dir, false; maxfreq=maxfreq_fft)

		if mol == 1 && load_params.how_many_molecules > 1
			plot!(p_tmp,
				xmirror=true,
				bottom_margin=2Plots.mm,
				top_margin=2Plots.mm,
				right_margin=8Plots.mm,
				left_margin=8Plots.mm,
			)
		elseif mol != load_params.how_many_molecules
			plot!(p_tmp,
				xlabel="",
				xticks=false,
				left_margin=8Plots.mm,
				right_margin=8Plots.mm,
				bottom_margin=2Plots.mm,
				top_margin=2Plots.mm
			)
			plot!(
				plot_mj_fft_tmp,
				xlabel="",
				xmirror=false,
				top_margin=8Plots.mm,
			)
			plot!(
				plot_ms_fft_tmp,
				xlabel="",
				xmirror=false,
				top_margin=8Plots.mm,
			)
		else
			plot!(p_tmp,
				bottom_margin=2Plots.mm,
				top_margin=2Plots.mm,
				right_margin=8Plots.mm,
				left_margin=8Plots.mm,
			)
			plot!(
				plot_mj_fft_tmp,
				xmirror=false,
				bottom_margin=16Plots.mm, # <--- TODO Change to 16Plots.mm if single molecule plotting
			)
			plot!(
				plot_ms_fft_tmp,
				xmirror=false,
				bottom_margin=16Plots.mm,
			)
		end
		push!(plot_vec_tmp, p_tmp)
		push!(fft_mj_vec_tmp, plot_mj_fft_tmp)
		push!(fft_ms_vec_tmp, plot_ms_fft_tmp)
	end
	# plot_final = plot(
	# 	plot_vec_tmp...,
	# 	layout=(load_params.how_many_molecules, 1),
	# 	legend_background_color="rgba(255, 255, 255, 0.6)",
	# 	left_margin=25Plots.mm,
	# 	right_margin=15Plots.mm,
	# 	size=(1400, 600*load_params.how_many_molecules),
	# )
	plot_final = plot(
		stack((fft_mj_vec_tmp, plot_vec_tmp, fft_ms_vec_tmp), dims=1)...,
		layout=grid(load_params.how_many_molecules, 3, widths=(1/6, 4/6, 1/6)),
		legend_background_color="rgba(255, 255, 255, 0.6)",
		# left_margin=15Plots.mm,
		# right_margin=35Plots.mm, 
		size=(1800, 1000*(load_params.how_many_molecules)),
	)

	savefig(joinpath(plot_png_dir, "03_quench_observables_m_ms_all.png"))
	savefig(joinpath(plot_pdf_dir, "03_quench_observables_m_ms_all.pdf"))

	plot_final
end

begin
	p_vec_tmp = []
	fft_vec_tmp = []
	for mol in 1:load_params.how_many_molecules
		color_1 = :cyan3
		color_2 = :magenta4
		lw = 2

		m_label = "m$(mol)"
		ms_label = "ms$(mol)"

		p_tmp = plot(
			dfs_observables[m_label].time_harmonic,
			dfs_observables[m_label].value + dfs_observables[ms_label].value,
			label=L"m_{j}^{(%$(mol))} + m_{s}^{(%$(mol))}",
			# ylabel=L"\langle m_{tot}^{(%$(mol))} \rangle =  \langle m_{j}^{(%$(mol))} + m_{s}^{(%$(mol))} \,\rangle",
			ylabel=L"\langle m_{tot}^{(%$(mol))} \rangle",
			xlabel="Time " * L"(2\pi/\omega_{\mathrm{trap}})",
			lw=lw,
			color=color_2,
			# y_guidefontcolor=color_1,
			# y_foreground_color_axis=color_1,
			# y_foreground_color_text=color_1,
			# y_foreground_color_border=color_1,
		)

		plot!(
			legend=false,
			xlim=xlim,
			frame=:box,
		)

		# savefig(joinpath(plot_png_dir, "quench_observables_mtot$(mol).png"))
		# savefig(joinpath(plot_pdf_dir, "quench_observables_mtot$(mol).pdf"))

		if mol != load_params.how_many_molecules
			plot!(p_tmp,
				title=L"\mathrm{(%$(lab_var))}~R_{ij} = %$(Int(load_params.geom_values_flat[1]))\,\mathrm{nm},\,\theta = %$(angle_var)^\circ",
				xmirror=true,
				bottom_margin=2Plots.mm,
				top_margin=20Plots.mm,
				right_margin=10Plots.mm,
				)
		else
			plot!(p_tmp,
			bottom_margin=15Plots.mm,
			top_margin=2Plots.mm,
			right_margin=10Plots.mm,
			xlabel="Time " * L"(2\pi/\omega_{\mathrm{trap}})",
		)
		end

		fft_plot_tmp = plot_fft_for_observable("mtot$(mol)", load_params.harmonic_frequency, dfs_observables, plot_dir, false; maxfreq=maxfreq_fft)

		if !(mol == 1 || mol == load_params.how_many_molecules)
			plot!(
				fft_plot_tmp,
				ylabel="",
				yticks=false,
				bottom_margin=2Plots.mm,
				top_margin=2Plots.mm,
				left_margin=10Plots.mm,
				right_margin=2Plots.mm,
			)
		elseif mol == load_params.how_many_molecules
			plot!(
				fft_plot_tmp,
				xmirror=false,
				)
		end
		fft_vec_tmp = push!(fft_vec_tmp, fft_plot_tmp)
		
		p_vec_tmp = push!(p_vec_tmp, p_tmp)
	end

	plot_final = plot(
		stack((p_vec_tmp, fft_vec_tmp), dims=1)...,
		layout=grid(load_params.how_many_molecules, 2, widths=(4/5, 1/5)),
		legend_background_color="rgba(255, 255, 255, 0.6)",
		# left_margin=15Plots.mm,
		# right_margin=35Plots.mm, 
		size=(1800, 300*(load_params.how_many_molecules)),
	)

	savefig(joinpath(plot_png_dir, "02_quench_observables_mtot_all.png"))
	savefig(joinpath(plot_pdf_dir, "02_quench_observables_mtot_all.pdf"))
end

begin 
	p_mtot_and_M_Ms = plot(
		dfs_observables["M"].time_harmonic,
		dfs_observables["M"].value,
		label=obs_rename_dict["M"],
		palette=:tab10,
	)

	plot!(
		p_mtot_and_M_Ms,
		dfs_observables["Ms"].time_harmonic,
		dfs_observables["Ms"].value,
		label=obs_rename_dict["Ms"],
		palette=:tab10,
	)

	for mol in 1:load_params.how_many_molecules
		color_1 = :cyan3
		color_2 = :magenta4
		lw = 2

		m_label = "m$(mol)"
		ms_label = "ms$(mol)"

		plot!(
			p_mtot_and_M_Ms,
			dfs_observables[m_label].time_harmonic,
			dfs_observables[m_label].value + dfs_observables[ms_label].value,
			label=L"m_{tot}^{(%$(mol))}",
			lw=lw,
			palette=:tab10,
		)

		plot_fft_for_observable("mtot$(mol)", load_params.harmonic_frequency, dfs_observables, plot_dir, true)
	end

	plot!(
		p_mtot_and_M_Ms,
		# title="Total Angular Momentum Projections vs Time",
		title=L"\mathrm{(%$(lab_var))}~R_{ij} = %$(Int(load_params.geom_values_flat[1]))\,\mathrm{nm},\,\theta = %$(angle_var)^\circ",
		xlabel="Time " * L"(2\pi/\omega_{\mathrm{trap}})",
		ylabel=L"\langle \Psi(t) | \hat{O}~|\Psi(t) \rangle",
		legend=:topright,
		xlim=xlim,
		legend_background_color="rgba(255, 255, 255, 0.6)",
		bottom_margin=15Plots.mm,
		size=(1600, 400),
	)

	savefig(joinpath(plot_png_dir, "01_quench_observables_all_mtot_and_M_Ms.png"))
	savefig(joinpath(plot_pdf_dir, "01_quench_observables_all_mtot_and_M_Ms.pdf"))

	p_mtot_and_M_Ms
end