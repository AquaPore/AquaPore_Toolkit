# =============================================================
#		module: plotHypix
# =============================================================
module plotHypix
	import  ..cst, ..kunsat, ..rootWaterUptake, ..tool, ..wrc, ..ΨminΨmax
	import Dates: value, DateTime


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT_HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_HYPIX(∑T_Reduced, ∑ΔQ_Obs_Reduced, clim, Date_Reduced, dateHypix, discret, hydro, hydroHorizon, i∑T_CalibrStart_Day, iMultistep, iScenario, N_iRoot, N_Layer, Nit_Reduced, Nz, optionHypix, paramHypix, pathInputHypix, pathOutputHypix, SiteName, veg, Z, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPrGross_Reduced, ΔPrThroughfall_Reduced, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRootDensity, ΔRunoff_Reduced, ΔSink_Reduced, θ_Reduced; obsθ=0, θobs_Reduced=[0;0], θsim_Aver=[0;0])

			if optionHypix.Plot_Hypix
				if optionHypix.θobs
					plotHypix.makkie.TIMESERIES(∑T_Reduced, ∑ΔQ_Obs_Reduced, Date_Reduced, dateHypix, discret, iMultistep, iScenario, Nit_Reduced, Nz, optionHypix, paramHypix, pathOutputHypix, pathInputHypix, SiteName, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPrGross_Reduced, ΔPrThroughfall_Reduced, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔSink_Reduced, θ_Reduced, θsim_Aver; obsθ=obsθ, θobs_Reduced=θobs_Reduced)
				else
					plotHypix.makkie.TIMESERIES(∑T_Reduced, ∑ΔQ_Obs_Reduced, Date_Reduced, dateHypix, discret, iMultistep, iScenario, Nit_Reduced, Nz, optionHypix, paramHypix, pathOutputHypix, pathInputHypix, SiteName, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPrGross_Reduced, ΔPrThroughfall_Reduced, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔSink_Reduced, θ_Reduced, θsim_Aver)
				end
			end

			if optionHypix.Plot_θprofile
				plotHypix.makkie.θPROFILE(∑T_Reduced, discret, iScenario, Nz, obsθ, optionHypix, paramHypix, pathOutputHypix, SiteName, θ_Reduced)
			end  # if: optionHypix.Plot_

			if optionHypix.Plot_θΨK
				plotHypix.θΨK(hydroHorizon, N_Layer, iMultistep, pathOutputHypix)
			end
			
			if optionHypix.Plot_Vegetation && optionHypix.RootWaterUptake
				plotHypix.VEG_FUNCTIONS(discret, iMultistep, N_iRoot, veg, Z, ΔRootDensity, pathOutputHypix)
			end
			
			if optionHypix.Plot_Interception
				plotHypix.plots.RAINFALL_INTERCEPTION(clim, i∑T_CalibrStart_Day, iMultistep, pathOutputHypix)
			end

			if optionHypix.Plot_HeatMap_YearMonth
				PathInput = pathOutputHypix.Table_Statistic * "_SeZaver_YearMonth_" * string(iScenario) * ".csv"
				PathOutput = pathOutputHypix.Plot_Heatmap_YearMonth * ".svg"

				plotHypix.makkie.HEATMAP(SiteName[iScenario], paramHypix, PathInput, PathOutput)
			end

			if optionHypix.Plot_HeatMap_MonthCombine
				PathInput = pathOutputHypix.Table_Statistic * "_SeZaver_YearMonth_Combine_" * string(iScenario) * ".csv"
				PathOutput = pathOutputHypix.Plot_Heatmap_MonthCombine * ".svg"

				plotHypix.makkie.HEATMAP(SiteName[iScenario], paramHypix, PathInput, PathOutput)
			end
			
			if  optionHypix.Plot_Sorptivity
				plotHypix.plots.PLOT_SORPTIVITY(hydro, iMultistep, optionHypix, pathOutputHypix)
			end
		return nothing
		end  # function: PLOT_HYPIX
	# ------------------------------------------------------------------



	# =============================================================
	#		module: makkie
	# =============================================================
	module makkie
		using CairoMakie
		using Dates, Statistics
		using DataFrames, CSV
		export θPROFILE, TIMESERIES

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θPROFILE(∑T_Reduced, discret, obsθ, optionHypix, paramHypix, θ_Reduced)
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θPROFILE(∑T_Reduced, discret, iScenario, Nz, obsθ, optionHypix, paramHypix, pathOutputHypix, SiteName, θ_Reduced)
			# PATH
				Path = pathOutputHypix.Plot_θprofile * ".svg"
				println("			 ~ ", Path, "~")
				rm(Path, force=true, recursive=true)

			# DEPTHS TO PLOT
				Zprofile = fill(0.0::Float64, Nz)
				for iZ=1:Nz
					Zprofile[iZ] = discret.Znode[iZ] / 100.0
				end

			# SELECTING PROFILE TO PLOT
				Nt = obsθ.Nit

			# INITIALIZING PLOT
				CairoMakie.activate!()
				Makie.inline!(true)

				Color_Hypix = [:red, :darkviolet, :orange,  :blue, :teal]

				Fig = Figure(resolution=(600,500))
				Title = SiteName[iScenario]

				Label_HyPix =fill("", Nt)
				Label_Hydrus =fill("", Nt)

				Ax1 = Axis(Fig[1,1], title=Title, xlabel= L"$\theta$  $[m^{3}  m^{-3}]$", ylabel= L"Z  $[cm]$", titlesize=25, xlabelsize=22, ylabelsize=22, xgridvisible=true, ygridvisible=false)

				# Ax2 = Axis(Fig[1,1],  font = "CMU Serif", titlesize=30, fontsize=16, xlabelsize=24, ylabelsize=24)

			# For every θprofile_Time
				for iT=1:Nt
					Tprofile = obsθ.∑T[iT]

					iTprofile = 1
	
					iTprofile = findfirst(x->x==Tprofile, ∑T_Reduced)
					
					if isnothing(iTprofile)
						println("Error θprofile_Time must be one of = $(obsθ.∑T)")
						error()
					end

					θprofile = θ_Reduced[iTprofile, 1:Nz]

					# PLOTTING
						# Label
							if paramHypix.ΔT_Output==3600.0 
								Label_HyPix[iT] = "HyP=" * string(ceil(Int, Tprofile / paramHypix.ΔT_Output)) * "Hour" 
							else
								Label_HyPix[iT] = "HyP_" * string(ceil(Int, Tprofile / paramHypix.ΔT_Output)) * "Day" 
							end

							if paramHypix.ΔT_Output==3600.0 
								Label_Hydrus[iT] = "HYD_" * string(ceil(Int, Tprofile / paramHypix.ΔT_Output)) * "Hour" 
							else
								Label_Hydrus[iT] = "HYD_" * string(ceil(Int, Tprofile / paramHypix.ΔT_Output)) * "Day" 
							end

					Plot2 = lines!(Ax1, obsθ.θobs[iT,1:Nz], -Zprofile, color=Color_Hypix[iT], linewidth=3, label=Label_Hydrus[iT])
					Plot1 = lines!(Ax1, θprofile, -Zprofile, color=Color_Hypix[iT], linewidth=2, label=Label_HyPix[iT], linestyle=:dash)
				end

				Leg = Legend(Fig[2,1], Ax1, framevisible=true, orientation=:horizontal, tellheight=true, nbanks=2, labelsize=14)
		
				trim!(Fig.layout)

				# axislegend()
				display(Fig)
				save(Path, Fig)

		return nothing
		end  # function: θPROFILE
		# ------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : TIMESERIES
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIMESERIES(∑T_Reduced::Vector{Float64}, ∑ΔQ_Obs_Reduced::Vector{Float64}, Date_Reduced::Vector{Dates.DateTime}, dateHypix, discret, iMultistep::Int64, iScenario::Int64, Nit_Reduced::Int64, Nz::Int64, optionHypix, paramHypix, pathOutputHypix, pathInputHypix, SiteName::Vector{Any}, ΔEvaporation_Reduced::Vector{Float64}, ΔPet_Reduced::Vector{Float64}, ΔPond_Reduced::Vector{Float64}, ΔPrGross_Reduced::Vector{Float64}, ΔPrThroughfall_Reduced::Vector{Float64}, ΔQ_Obs_Reduced, ΔQ_Reduced::Matrix{Float64}, ΔRunoff_Reduced::Vector{Float64}, ΔSink_Reduced::Vector{Float64}, θ_Reduced::Matrix{Float64}, θsim_Aver; obsθ=0, θobs_Reduced=[0 ; 0])

			# Dimensions of figure
				Height= 600
				Width = 1000

			# Number of days to print 
				Nticks = 24

			# PATH
				Path = pathOutputHypix.Plot_HypixTime * "_" * string(iMultistep) * ".svg"
				rm(Path, force=true, recursive=true)

			# STYLE
				Style_Color = [:red, :darkviolet, :orange, :teal, :blue,:brown]

			# TICKS
				Date_SimStart   = dateHypix.Date_SimStart  # since we need to compute the culmulative of the 1rst day
				Date_End_Calibr = dateHypix.Date_End
				if optionHypix.θobs
					Date_SimStart   = max(obsθ.Date[1], dateHypix.Date_SimStart)  # since we need to compute the culmulative of the 1rst day
					Date_End_Calibr = min(obsθ.Date[end], dateHypix.Date_End)
				else
					Date_SimStart   = dateHypix.Date_SimStart  # since we need to compute the culmulative of the 1rst day
					Date_End_Calibr = dateHypix.Date_End
				end
				
			# DAYS PLOT
				ΔDays = floor(Int, Dates.value((Date_End_Calibr - Date_SimStart)) / ((Nticks - 1) * paramHypix.ΔT_Output * 1000))::Int64

				# Putting a monthly dates
					Ndates = length(∑T_Reduced)

					Ndates_Reduced = floor(Int, ∑T_Reduced[Ndates] / (ΔDays * paramHypix.ΔT_Output)) + 1
					
					Date_Reduced2 = fill("", Ndates_Reduced+1)
					∑T_Reduced2 = fill(0::Int64, Ndates_Reduced+1)

					# Reducing dates
					∑T_Reduced2[1] = ∑T_Reduced[1]	
	
					Date_Reduced2[1]= Dates.format(Date_Reduced[1], "d u Y")
					iGood = 2
					for i=1:Ndates
						if ∑T_Reduced[i] -  ∑T_Reduced[1] ≥ ((iGood - 1) * ΔDays * paramHypix.ΔT_Output)

							∑T_Reduced2[iGood] = ∑T_Reduced[i]	
	
							Date_Reduced2[iGood]= Dates.format(Date_Reduced[i], "d u Y")

							iGood += 1
						end # if
					end # for
			
			# PLOTTING
			# , resolution = (3000, 2500)
				CairoMakie.activate!(type = "svg")
				# fontsize=40,
				Fig = Figure(font="Sans", titlesize=50,  xlabelsize=30, ylabelsize=30, labelsize=30, fontsize=30)
			
			# Plot Climate -------------------------------------------------	
			iSubplot = 0
			if optionHypix.Plot_Climate
				iSubplot += 1
				Title = "©HyPix  " * SiteName[iScenario]

				Axis1 = Axis(Fig[iSubplot,1], title=Title, ylabel= L"$\Delta Fluxes$ $[mm$ $day ^{-1}]$", xgridvisible=true, ygridvisible=false, height=Height, width=Width, yticks = LinearTicks(10))

					hidexdecorations!(Axis1, grid=false, ticks=true, ticklabels=true)

					Axis1.xticks = (∑T_Reduced2[1:iGood], string.(Date_Reduced2[1:iGood]))

					xlims!(Axis1, ∑T_Reduced[1], ∑T_Reduced[Nit_Reduced]) 
					# ylims!(Axis1)
	
					Label3= L" $\Delta PrTop$"
					Plot_Climate3 = barplot!(Axis1, ∑T_Reduced[1:Nit_Reduced], ΔPrGross_Reduced[1:Nit_Reduced], strokecolor=:blue, strokewidth=1.5, color=:blue)

					Label5=L"$\Delta PrThrough"	
					Plot_Climate5 = barplot!(Axis1, ∑T_Reduced[1:Nit_Reduced], ΔPrThroughfall_Reduced[1:Nit_Reduced], strokecolor=:chartreuse3, strokewidth=1.5, color=:olive)

					Label1=L"$\Delta Runoff$"	
					Plot_Climate1 = barplot!(Axis1, ∑T_Reduced[1:Nit_Reduced], ΔRunoff_Reduced[1:Nit_Reduced], strokewidth=1.5, strokecolor=:magenta2,  color=:magenta2)

					Label4=L"$\Delta PrSoil$"					
						Plot_Climate4 = barplot!(Axis1, ∑T_Reduced[1:Nit_Reduced], ΔQ_Reduced[1:Nit_Reduced, 1], strokecolor=:cyan, strokewidth=1.5, color=:cyan)

					Label2=L"$Hpond$"
						Plot_Climate2 = barplot!(Axis1, ∑T_Reduced[1:Nit_Reduced], -ΔPond_Reduced[1:Nit_Reduced], strokewidth=0.5, color=:chocolate1, strokecolor=:chocolate1)
					
					Legend(Fig[iSubplot,2], [Plot_Climate3, Plot_Climate5, Plot_Climate4, Plot_Climate2, Plot_Climate1], [Label3, Label5, Label4, Label2, Label1], framevisible=true, tellwidth=true, tellheight=true, margin=(10, 10, 10, 10))

					
				# GROUND WATER RECHARGE AXIS 3 ===================================================================										
					# trim!(Fig.layout)

					# Label7= L"$\Delta Q$"
					iSubplot += 1
					Axis3  = Axis(Fig[iSubplot, 1], xgridvisible=true, ygridvisible=false, height=Height*0.6, width=Width, ylabel= L"$\Delta Q$ $[mm$ $day ^{-1}]$")

					if !(isempty(pathInputHypix.Drainage[iScenario]))
						Axis3b = Axis(Fig[iSubplot, 1], yaxisposition=:right, xgridvisible=true, ygridvisible=false, ylabel= L"$\sum \Delta Q$ $[mm$ $day ^{-1}]$")
						hidespines!(Axis3b)
						hidexdecorations!(Axis3b)
					end

					Axis3.xticks = (∑T_Reduced2[1:iGood], string.(Date_Reduced2[1:iGood]))
						xlims!(Axis3, ∑T_Reduced[1],∑T_Reduced[Nit_Reduced])
						# ylims!(Axis3, -maximum(ΔQ_Reduced[1:Nit_Reduced, Nz+1])-eps(), 0)

						hidexdecorations!(Axis3, grid=false, ticks=true, ticklabels=true)
						Label_7a = L"$\Delta Q hypix$"
						Plot_Climate7a = barplot!(Axis3, ∑T_Reduced[1:Nit_Reduced], -ΔQ_Reduced[1:Nit_Reduced, Nz+1], strokecolor=:violetred, strokewidth=1, color=:violet)

						if !(isempty(pathInputHypix.Drainage[iScenario]))
							Label_7b = L"$\Delta Q _{Obs}$"
							Plot_Climate7b = barplot!(Axis3, ∑T_Reduced[1:Nit_Reduced], -ΔQ_Obs_Reduced[1:Nit_Reduced], strokecolor=:violet, strokewidth=1, color=:red1, label=Label_7b)

							Label_7c = L"$ \sum \Delta Q _{Obs}$"
							Plot_Climate7c = lines!(Axis3b, ∑T_Reduced[1:Nit_Reduced], -∑ΔQ_Obs_Reduced[1:Nit_Reduced], color=:firebrick, linewidth=2.5, label=Label_7c)

							# To start from 0
							∑ΔQ_Reduced = fill(0.0::Float64, Nit_Reduced)
							for iT=2:Nit_Reduced
								∑ΔQ_Reduced[iT] = ∑ΔQ_Reduced[iT-1] + ΔQ_Reduced[iT, Nz+1]  
							end
							Label_7d = L"$\sum \Delta Q _{Hypix}$"
							Plot_Climate7d = lines!(Axis3b, ∑T_Reduced[1:Nit_Reduced], -∑ΔQ_Reduced[1:Nit_Reduced], color=:orchid1, label=Label_7d, linestyle=:dash, linewidth=2.5)

							Legend(Fig[iSubplot,2], [Plot_Climate7a, Plot_Climate7b, Plot_Climate7c, Plot_Climate7d], [Label_7a, Label_7b, Label_7c, Label_7d], framevisible=true, tellwidth=true, tellheight=true, margin=(10, 10, 10, 10))
						else
							Legend(Fig[iSubplot,2], [Plot_Climate7a], [Label_7a], framevisible=true, tellwidth=true, tellheight=true, margin=(10, 10, 10, 10))
						end



						# Leg = Legend(Fig[iSubplot:iSubplot-1,2], Axis3, framevisible=true, tellwidth=true, tellheight=true, orientation = :horizontal)
						# Leg = Legend(Fig[iSubplot:iSubplot+1,2], Axis3b, framevisible=true, tellwidth=true, tellheight=true, orientation = :horizontal)
					
						# Legend(Fig[iSubplot,2], [Plot_Climate7], [Label7], framevisible=true, tellwidth=true, tellheight=true, margin = (10, 10, 10, 10))

						iSubplot += 1
						# trim!(Fig.layout)
					

					# Plot_Climate = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], (ΔSink_Reduced[1:Nit_Reduced].-ΔEvaporation_Reduced[1:Nit_Reduced]), colour=:blue, label=L"$\Delta Rwu$")
			end # if: optionHypix.Plot_Climate


			# PLOT θ Axis 4 ===================================================================		
			if optionHypix.Plot_θ		
				iSubplot += 1
				Axis4 = Axis(Fig[iSubplot,1], ylabel=L"$\theta$ $[mm^3 mm^{-3}]$", xgridvisible=true, ygridvisible=false, height=Height, width=Width)

				Axis4.xticks = (∑T_Reduced2[1:iGood], string.(Date_Reduced2[1:iGood]))
				Axis4.xticklabelrotation = π/3
				hidexdecorations!(Axis4, grid=false, ticks=true, ticklabels=true)

				xlims!(Axis4, ∑T_Reduced[1],∑T_Reduced[Nit_Reduced])

				if optionHypix.θobs
					Nbanks = 2

					# Observation θplot obs
					for iZobs = 1:obsθ.Ndepth
						# LABEL
							Label_Obs = "θobs_" * string(Int(floor(obsθ.Z[iZobs]))) * "mm"
							Label_Sim = "θhypix_" * string(Int(floor((discret.Znode[obsθ.ithetaObs[iZobs]])))) * "mm"

							if optionHypix.θobs
								if optionHypix.θavr_RootZone
									Plot_θsim = lines!(Axis4, ∑T_Reduced[1:Nit_Reduced], θsim_Aver[1:Nit_Reduced,1], linewidth=1.5, color=:blue, label=L"$\theta _{Hypix}$", linestyle=:dash)

									# Correcting θbs
										θsim_Aver_Mean    = mean(θsim_Aver[1:Nit_Reduced,1])
										θobs_Reduced_Mean = mean(filter(!isnan, θobs_Reduced[1:Nit_Reduced, 1]))
										θobs_Corrected    = θobs_Reduced[1:Nit_Reduced, iZobs] .+ (θsim_Aver_Mean - θobs_Reduced_Mean)
										θobs_Corrected_Min = minimum(filter(!isnan, θobs_Corrected))
										if θobs_Corrected_Min < 0.0
											@. θobs_Corrected = θobs_Corrected + θobs_Corrected_Min
										end

									Axis4.xticks = (∑T_Reduced2[1:iGood], string.(Date_Reduced2[1:iGood]))
									hidexdecorations!(Axis4, grid=false, ticks=true, ticklabels=true)

									Plot_θobs = lines!(Axis4, ∑T_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, iZobs], linewidth=1.5, color=:red, label=L"$\theta _{Obs}$")

									Plot_θobs_Corrected = lines!(Axis4, ∑T_Reduced[1:Nit_Reduced], θobs_Corrected, linewidth=1.5, color=:violet, label=L"$\theta _{ObsCorected}$")

									Nbanks = 1
								else
									Axis4.xticks = (∑T_Reduced2[1:iGood], string.(Date_Reduced2[1:iGood]))
									Axis4.xticklabelrotation = π/2

									# Plot_θobs = lines!(Axis4, ∑T_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, iZobs], linewidth=1.5, color=Style_Color[iZobs], label=L"$\theta _{Hypix}$")
									Plot_θobs = lines!(Axis4, ∑T_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, iZobs], linewidth=1.5, color=Style_Color[iZobs], label=Label_Obs)

									Plot_θsim = lines!(Axis4, ∑T_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced, obsθ.ithetaObs[iZobs]], linewidth=1.5, color=Style_Color[iZobs], label=Label_Sim, linestyle=:dash)

									Nbanks = obsθ.Ndepth
								end
							end
					end # loop
				else
					# Depths to plot
					Nθplots = 5
					Nθplots = min(Nθplots, Nz)
					ΔCell = Int64(floor((Nz-1) / Nθplots))
					Nbanks = Nθplots-1

						CellPlot = [1]
						for i = 1:Nθplots-1
							∑ΔCell = min(CellPlot[end] + ΔCell, Nz)
							append!(CellPlot, ∑ΔCell)
						end 

					# Plotting
					for iZobs = 1:Nθplots
						# LABEL
							Label_Sim = "θhypix" * string(Int(floor((discret.Znode[CellPlot[iZobs]])))) * "mm"

						# PLOTS
							Plot_θsim = lines!(Axis4, ∑T_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced, CellPlot[iZobs]], linewidth=1.5, color=Style_Color[iZobs], label=Label_Sim, linestyle=:dash)
					end # loop						
				end # if optionHypix.θobs

				Fig[iSubplot, 2] = Legend(Fig, Axis4, framevisible=true, tellwidth=true, tellheight=true, margin = (10, 10, 10, 10))

			end # if: optionHypix.Plot_θ


			if optionHypix.θavr_RootZone && optionHypix.θobs
				iSubplot += 1

				Axis5 = Axis(Fig[iSubplot,1], ylabel=L"$\theta _{MeanAdj}$ $[mm^3 mm^{-3}]$", xgridvisible=true, ygridvisible=false, height=Height, width=Width, xlabelsize=30)

				xlims!(Axis5, ∑T_Reduced[1],∑T_Reduced[Nit_Reduced])


				if !optionHypix.Plot_Etp
					Axis5.xticks = (∑T_Reduced2[1:iGood], string.(Date_Reduced2[1:iGood]))
					Axis5.xticklabelrotation = π/3
				else
					hidexdecorations!(Axis5, grid=false, ticks=true, ticklabels=true)
				end
				
				iZobs = 1

				θsim_Aver_Mean    = mean(θsim_Aver[1:Nit_Reduced,1])
				θobs_Reduced_Mean = mean(filter(.!isnan, θobs_Reduced[1:Nit_Reduced, 1]))

				Plot_θobs = lines!(Axis5, ∑T_Reduced[1:Nit_Reduced], θobs_Reduced[1:Nit_Reduced, iZobs] .- θobs_Reduced_Mean, linewidth=3.0, color=:red, label= L"$\theta MeanAdjObs$")
				
				Plot_θsim = lines!(Axis5, ∑T_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced, obsθ.ithetaObs[iZobs]] .- θsim_Aver_Mean , linewidth=3.0, color=:blue, label=L"$\theta MeanAdjHypix$", linestyle = :dash)

				Plot_Line = lines!(Axis5, ∑T_Reduced[1:Nit_Reduced], 0.0 .* ∑T_Reduced[1:Nit_Reduced] , linewidth=2, color=:silver,  linestyle = :dashdot)

				Fig[iSubplot, 2] = Legend(Fig, Axis5, framevisible=true, tellwidth=true, tellheight=true, margin = (10, 10, 10, 10))
			end

		# POTENTIAL EVAPOTRANSPIRATION AXIS 2 ===================================================================		
			if optionHypix.Plot_Etp
				iSubplot += 1
				Axis2  = Axis(Fig[iSubplot, 1], yticklabelcolor=:black, yaxisposition=:left, rightspinecolor=:black, ytickcolor=:black, ylabel= L"$\Delta ET$ $[mm$ $day ^{-1}]$", xgridvisible=true, ygridvisible=false, height=Height, width=Width)

				Axis2.xticks = (∑T_Reduced2[1:iGood], string.(Date_Reduced2[1:iGood]))
				Axis2.xticklabelrotation = π/3

					xlims!(Axis2, ∑T_Reduced[1],∑T_Reduced[Nit_Reduced])
					ylims!(Axis2, low=0)

					Label4 =L"$PotEvapoTransp$"
					Plot_Climate4 = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], ΔPet_Reduced[1:Nit_Reduced], linewidth=2, colour=:darkgreen)

					Label5 = L"$\Delta EvapoTransp$"
					Plot_Climate5 = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], ΔSink_Reduced[1:Nit_Reduced], linewidth=2, colour=:red)

					Label6=L"$\Delta Evapo$"
					Plot_Climate6 = lines!(Axis2, ∑T_Reduced[1:Nit_Reduced], ΔEvaporation_Reduced[1:Nit_Reduced], linewidth=2, colour=:purple4)

					Legend(Fig[iSubplot,2], [Plot_Climate4, Plot_Climate5, Plot_Climate6], [Label4, Label5, Label6], framevisible=true, tellwidth=true, tellheight=true, margin = (10, 10, 10, 10))
				end

		colgap!(Fig.layout, 15)
		rowgap!(Fig.layout, 15)
		trim!(Fig.layout)
		
		resize_to_layout!(Fig)
		display(Fig)
		save(Path, Fig, pt_per_unit=0.5)
		println("			 ~ ", Path, "~")
		
		return nothing
		end  # function: TIMESERIES
		
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HEATMAP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HEATMAP(SiteName₁, paramHypix, PathInput, PathOutput)
			
			Df = DataFrame(CSV.File(PathInput))

			Df_θZₐᵥₑᵣ =  "Z" .* string.(paramHypix.table.θZₐᵥₑᵣ) .*"mm" .* "_mean"
		
			Header_θZₐᵥₑᵣ = reverse(string.(paramHypix.table.θZₐᵥₑᵣ))

			Nit = length(Df[!,Df_θZₐᵥₑᵣ[1]])
			N_θZₐᵥₑᵣ = length(Df_θZₐᵥₑᵣ)

			Width = max(Nit * 40, 800)
			Height= 500

		# PLOTTING  ====
			Xticks = string.(Df[!,:Year]) .* "_" .* Dates.monthabbr.(Df[!,:Month])
			Fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), font="CMU Serif", titlesize=50, fontsize=25, xlabelsize=11, ylabelsize=25, labelsize=25)
			
			Title = "©HyPix,    SiteName = $SiteName₁" 
			Ax1 = Axis(Fig[1,1], title=Title, titlesize=40, xticks =(1:Nit, Xticks), ylabel= L"$\Delta Flux$ [mm]", height=Height, width=Width, xgridvisible = true, ygridvisible = false, yticks = LinearTicks(10))
				xlims!(Ax1, 1, Nit)

				hidexdecorations!(Ax1, grid=false, ticks=true, ticklabels=true)
				
				barplot!(Ax1, 1:Nit, Df[!,:ΔPrThroughfall_sum], color=:blue, strokecolor=:black, strokewidth=1.5,  label= L"$\Delta PrThrough$ [mm]")

				barplot!(Ax1, 1:Nit, Df[!,:ΔRunoff_sum], color=:magenta2, strokecolor=:black, strokewidth=1.5,  label= L"$\Delta Runoff$ [mm]")
				
				barplot!(Ax1, 1:Nit, - Df[!,:ΔQ_sum], color=:red, strokecolor=:black, strokewidth=1.5,  label=L"$\Delta Q$ [mm]")

				lines!(Ax1, 1:Nit, Df[!,:ΔSink_sum], linewidth=2.5, colour=:olivedrab1, label= L"$\Delta EvapTransp$ [mm]")

				Leg = Legend(Fig[1,2], Ax1, framevisible=true, tellheight=true, tellwidth=true, labelsize=25)
			
			Yticks = Header_θZₐᵥₑᵣ
			Ax2 = Axis(Fig[3, 1], xticks = (1:Nit, Xticks), yticks = (1:N_θZₐᵥₑᵣ, Yticks), height=Height, width=Width, ylabel= L"$Z$ [mm]")
				Ax2.xticklabelrotation = π / 3

				θZₐᵥₑᵣ_2D = Matrix{Float64}(Df[!, Df_θZₐᵥₑᵣ])
				θZₐᵥₑᵣ_2D = reverse(θZₐᵥₑᵣ_2D, dims=2)

				Hmap = heatmap!(Ax2, θZₐᵥₑᵣ_2D, colorrange=(0.0, 1.0), colormap = Reverse(:linear_kbc_5_95_c73_n256))

				Colorbar(Fig[3, 2], Hmap; width=20, ticks = 0:0.25:1, label=L"Water Filled Pore Space:  $\theta$ \ $\theta _S$  [-]", labelsize=25)

			# Writting the text nand change colour of the text
			for iT in 1:Nit, iZ in 1:N_θZₐᵥₑᵣ
				txtcolor = θZₐᵥₑᵣ_2D[iT, iZ] > 0.5 ? :white : :black
				text!(Ax2, "$(round(θZₐᵥₑᵣ_2D[iT,iZ], digits = 2))", position = (iT, iZ),
					color = txtcolor, align = (:center, :center), fontsize=18)
			end

			Ax2.xticklabelalign = (:right, :center)

			colgap!(Fig.layout, 15)
			rowgap!(Fig.layout, 15)
			resize_to_layout!(Fig)
			trim!(Fig.layout)
			display(Fig)
			save(PathOutput, Fig, pt_per_unit=0.5) # size = 600 x 450 pt
		
	return nothing
	end  # function: HEATMAP
# ------------------------------------------------------------------	

	end  # module: makkie

end  # module plotHypix
# ............................................................


	# ========================================
	# PLOTTING HYDRAULIC RELATIONSHIP FOR EVERY HORIZON
	# ======================================== 
	# function θΨK(hydroHorizon, N_Layer, iMultistep, pathOutputHypix)

	# 	# Deriving the Min and Max Ψ from principals of soil physics
	# 	Ψ_Min_Horizon = fill(0.0::Float64, N_Layer)
	# 	Ψ_Max_Horizon = fill(0.0::Float64, N_Layer)
	# 	for iZ=1:N_Layer
	# 		Ψ_Max_Horizon[iZ], Ψ_Min_Horizon[iZ] = ΨminΨmax.ΨMINΨMAX(hydroHorizon.θs[iZ], hydroHorizon.θsMacMat[iZ], hydroHorizon.σ[iZ], hydroHorizon.σmac[iZ], hydroHorizon.Ψm[iZ], hydroHorizon.ΨmMac[iZ])
	# 	end  # for iZ=1:N_Layer
		
	# 	# PREPARING THE DATA
	# 		N_Se = 1000
	# 		local Ψplot = exp.(range(log(minimum(Ψ_Min_Horizon[1:N_Layer])), stop = log(maximum(Ψ_Max_Horizon[1:N_Layer])), length=N_Se)) 

	# 		local θplot    = fill(0.0::Float64, N_Se)
	# 		local Kplot    = fill(0.0::Float64, N_Se)
	# 		local ∂θ∂Ψplot = fill(0.0::Float64, N_Se)
	# 		local ∂K∂Ψplot = fill(0.0::Float64, N_Se)

	# 		Plot_θΨK = PGFPlots.GroupPlot(4, 100, groupStyle = "horizontal sep = 3.5cm, vertical sep = 3.5cm")

	# 	# FOR EVERY HORIZON
	# 	for iZ = 1:N_Layer
			
	# 		for iΨ = 1:N_Se
	# 			if Ψ_Max_Horizon[iZ] ≥ Ψplot[iΨ] ≥ Ψ_Min_Horizon[iZ]
	# 				θplot[iΨ]    = wrc.Ψ_2_θ(optionₘ,Ψplot[iΨ], iZ, hydroHorizon)
					
	# 				Kplot[iΨ]    = kunsat.KUNSAT_θΨSe(optionₘ, Ψplot[iΨ], iZ, hydroHorizon)
					
	# 				∂θ∂Ψplot[iΨ] = wrc.∂θ∂Ψ(optionₘ, Ψplot[iΨ], iZ, hydroHorizon)

	# 				∂K∂Ψplot[iΨ] = kunsat.∂K∂ΨMODEL(optionₘ, Ψplot[iΨ], iZ, hydroHorizon)
	# 			else
	# 				θplot[iΨ]    = NaN
					
	# 				Kplot[iΨ]    = NaN
					
	# 				∂θ∂Ψplot[iΨ] = NaN

	# 				∂K∂Ψplot[iΨ] = NaN
	# 			end
	# 		end # for iΨ

	# 		Θs_Max = maximum(hydroHorizon.θs[1:N_Layer]) + 0.05
	# 		Ks_Min = 10.0 ^ -7 * cst.MmS_2_CmH
	# 		Ks_Max = maximum(hydroHorizon.Ks[1:N_Layer]) * cst.MmS_2_CmH * 1.1

	# 		Title =" $(pathOutputHypix.SiteName_Hypix)  Layer = $(iZ)"
		
	# 	# Plot 1: θΨ
	# 		Plot_θΨ = PGFPlots.Plots.Linear(log.(Ψplot) , θplot, style=" smooth, blue, very thick", mark="none", legendentry=L"$ \theta ( \Psi ) $")

	# 		Plot_hydro = [Plot_θΨ]

	# 		push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xlabel=L"$ Ln \ \Psi [mm]$", ylabel=L"$ \theta \ [mm{^3} \ mm^{-3}]$", xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), ymin=0.0, ymax=Θs_Max, legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

	# 	# Plot 2: Kplot(Ψplot)
	# 		Plot_Kθ = PGFPlots.Plots.Linear(log.(Ψplot), Kplot .* cst.MmS_2_CmH, style=" smooth, red, very thick", mark="none", legendentry=L"$ K_{unsat} \ ( \Psi ) $")

	# 		Plot_hydro = [Plot_Kθ]
			
	# 		push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title,  xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), ymax=Ks_Max, ymode="log", xlabel=L"$Ln \  \Psi [mm]$", ylabel=L"$ K_{unsat} \ [cm \ h^{-1}]$", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

	# 	# Plot 3: ∂θ∂Ψplot
	# 		Plot_∂θ∂Ψ = PGFPlots.Plots.Linear(log.(Ψplot), ∂θ∂Ψplot , style=" smooth, green, very thick", mark="none", legendentry=L"$ \partial \theta \partial \Psi $")

	# 		Plot_hydro = [Plot_∂θ∂Ψ]

	# 		push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), xlabel=L"$Ln \  \Psi [mm] $", ylabel=L"$ \partial \theta \partial \Psi $", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

	# 	# Plot 4: ∂K∂Ψplot
	# 		Plot_∂K∂Ψ = PGFPlots.Plots.Linear(log.(Ψplot), ∂K∂Ψplot, style=" smooth, teal, very thick", mark="none", legendentry=L"$ \partial K \partial \Psi $")

	# 		Plot_hydro = [Plot_∂K∂Ψ]

	# 		push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), xlabel=L"$Ln \  \Psi \ [mm]$", ylabel=L"$ \partial K \partial \Psi $", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

	# 	end #iZ ............................................

	# 	Path = pathOutputHypix.Plot_Hypix_θΨK * "_" * string(iMultistep) * ".svg"
	# 	PGFPlots.save(Path, Plot_θΨK) 
	# end # function θΨK


		# ............................................................

			# =============================================================
			#		module: plots
			# =============================================================
			# module plots
			# import ...sorptivity, ..wrc, ..cst
			# export PLOT_SORPTIVITY

			# 	using Plots.PlotMeasures, LaTeXStrings
			# 	using Plots
			# 	using Dates
				
			# 	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# 	#		FUNCTION : PLOT_SORPTIVITY
			# 	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# 	function PLOT_SORPTIVITY(hydro, iMultistep, optionHypix, pathOutputHypix)
			# 		println("  ==  START: PLOT_SORPTIVITY_SeIni  ==")

			# 		# Setting the range of values for Se
         #          Se_Ini         = collect(0.0:0.001:1.0)
         #          N_SeIni        = length(Se_Ini)
         #          Sorptivity_Mod = fill(0.0::Float64, (N_SeIni))
         #          θini          = fill(0.0::Float64, (N_SeIni))

			# 		for iSeIni=1:N_SeIni
			# 			θini[iSeIni] = wrc.Se_2_θ(Se_Ini[iSeIni], 1, hydro)

			# 			Sorptivity_Mod[iSeIni] = sorptivity.SORPTIVITY(θini[iSeIni], 1, hydro, optionHypix) 
			# 		end
					
			# 		# PLOTTING ====================	
			# 			Plot1=Plots.plot(layout=1)

			# 			Title =" $(pathOutputHypix.SiteName_Hypix)"

			# 			Plots.plot!(Plot1, Se_Ini[1:N_SeIni] , Sorptivity_Mod[1:N_SeIni], framestyle = [:box :semi :origin :zerolines :grid :true], xlabel=L"Initial \ Se \ [-]", ylabel=L"Sorptivity \  [ \ mm \ \sqrt s \ ]", label="", grid=false) 
					
			# 			Path =pathOutputHypix.Plot_Sorptivity  * "_" * string(iMultistep) * ".svg"

			# 			Plots.savefig(Plot1, Path)

			# 			println("			 ~ ", Path, "~")

			# 	end  # function: PLOT_SORPTIVITY
