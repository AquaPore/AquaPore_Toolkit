# =============================================================
#		MODULE: plot
#
# =============================================================
module plot
# =============================================================
#		MODULE: lab
# =============================================================
module lab
	import ...cst, ...kunsat, ...wrc, ...θψ_2_KsψModel
	using CairoMakie, ColorSchemes

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYDROPARAM
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYDROPARAM(hydro, hydroOther, IdSelect, K_KΨobs, NiZ, N_KΨobs, N_θΨobs, optim, option, param, path, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs; N_Se=1000)
				println("  ==  START: Plotting HydroParam  ==")

				# ===================== DATA =====================
            θ_Sim      = fill(0.0,N_Se)
            Kunsat_Sim = fill(0.0,N_Se)
            KsMat      = fill(0.0, param.globalparam.N_iZ_Plot_End)

				Ψ_θΨobs_Min = 0.0
				for iZ = param.globalparam.N_iZ_Plot_Start:param.globalparam.N_iZ_Plot_End

					Ψ_θΨobs_Max = maximum(Ψ_θΨobs[iZ,N_θΨobs[iZ]]) + 100000.0

					Ψ_Sim = expm1.(range(log1p(Ψ_θΨobs_Min), stop=log1p(Ψ_θΨobs_Max), length=N_Se)) 

					KsMat[iZ] =hydro.Ks[iZ] * min(max((hydro.θsMacMat[iZ] - hydro.θr[iZ]) / (hydro.θs[iZ] - hydro.θr[iZ]), 0.0), 1.0)

					θ_θΨobs_Max = hydro.Φ[iZ]

					# Simulated 
						for iΨ = 1:N_Se
							θ_Sim[iΨ] = wrc.Ψ_2_θ(option.hydro, Ψ_Sim[iΨ], iZ, hydro)
							Kunsat_Sim[iΨ] = kunsat.KUNSAT_θΨSe(option.hydro, Ψ_Sim[iΨ], iZ, hydro)
						end # iΨ = 1:N_Se

					# _______________________ START: Plotting _______________________
								
					Fig = Figure(resolution = (2500, 1000),  font="Sans", fontsize=16)

					Title = "iZ= $(IdSelect[iZ]) " * "θ(Ψ) Nse_θΨ=" * string(round(hydroOther.Nse_θΨ[iZ], digits=2)) * "; Nse_KΨ=" * string(round(hydroOther.Nse_KΨ[iZ], digits=2)) * "; Wilmot_θΨ=" *  string(round(hydroOther.NseWilmot_θΨ[iZ],digits=2)) * "; Wilmot_KΨ=" * string(round(hydroOther.NseWilmot_KΨ[iZ], digits=2))
					
					#  == Plot_θ_Ψ  ==
						Axis1 = Axis(Fig[1,1], title=Title, titlesize=24, xlabel="ln(1 + Ψ) [kPa]", ylabel="θ [mm³ mm⁻³]", xlabelsize=10, backgroundcolor=:white)

						xlims!(Axis1, log1p.(cst.Mm_2_kPa * Ψ_θΨobs_Min), log1p.(cst.Mm_2_kPa * Ψ_θΨobs_Max ))
						ylims!( Axis1, 0.0, max(hydro.Φ[iZ], maximum(θ_θΨobs[iZ,1:N_θΨobs[iZ]])) )

						Axis1.xticks = (log1p.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]]), string.( floor.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]], digits=1)))

						Fig_θΨobs = scatter!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_θΨobs[iZ,1:N_θΨobs[iZ]]), Float64.(θ_θΨobs[iZ,1:N_θΨobs[iZ]]), color=:red, markersize=25, marker = '■')

						Fig_θΨsim = lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_Sim[1:N_Se], color=:blue, linewidth=3)

						lines!(Fig[1,1], [Point(log1p(cst.Mm_2_kPa * hydro.ΨmacMat[iZ]), 0), Point(log1p(cst.Mm_2_kPa * hydro.ΨmacMat[iZ]), hydro.θsMacMat[iZ])], color=:brown, linewidth=3)

						lines!(Fig[1,1], [Point(log1p(0.0), hydro.θsMacMat[iZ]), Point(log1p(cst.Mm_2_kPa * hydro.ΨmacMat[iZ]), hydro.θsMacMat[iZ])], color=:brown, linewidth=3)
				
						Fig_TotalPorosity = scatter!(Fig[1,1], [log1p.(cst.Mm_2_kPa .* 0.0)], [hydro.Φ[iZ]], color=:green, markersize=25, marker ='●')

					# == Plot_K_Ψ  ==
					# If Ks is not computed it is computed from Ks_Model

						Axis2 = Axis(Fig[1,2], title="K(Ψ)", titlesize=24, xlabel = "ln(1 + Ψ) [kPa]", ylabel = "ln (1 + K (Ψ)) [mm h⁻¹]")

						xlims!(Axis2, log1p.(cst.Mm_2_kPa * Ψ_θΨobs_Min), log1p.(cst.Mm_2_kPa * Ψ_θΨobs_Max))

						ylims!(Axis2, 0.0, log1p(hydro.Ks[iZ]*cst.MmS_2_MmH))

						Axis2.xticks = (log1p.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]]), string.(floor.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]], digits=1)))
						Yticks = 1:1:6
						Axis2.yticks = (Yticks,string.(Yticks))

						if option.data.Kθ
							Fig_Kθobs = scatter!(Fig[1,2], log1p.(Ψ_KΨobs[iZ,1:N_KΨobs[iZ]].*cst.Mm_2_kPa), log1p.(K_KΨobs[iZ,1:N_KΨobs[iZ]].*cst.MmS_2_MmH), color=:red, markersize=25, marker = '■')
						end

						Fig_Kθsim = lines!(Fig[1,2], log1p.(Ψ_Sim[1:N_Se].*cst.Mm_2_kPa), log1p.(Kunsat_Sim[1:N_Se] .* cst.MmS_2_MmH), color=:blue, linewidth=3)

						Fig_Ks = scatter!(Fig[1,2], [log1p.(cst.Mm_2_kPa .* 0.0)], [log1p(hydro.Ks[iZ] * cst.MmS_2_MmH)], color=:green, markersize=25, marker ='●')

						lines!(Fig[1,2], [ Point(log1p(cst.Mm_2_kPa * hydro.ΨmacMat[iZ]), 0) , Point(log1p(cst.Mm_2_kPa * hydro.ΨmacMat[iZ]), log1p(KsMat[iZ]* cst.MmS_2_MmH))], color=:brown, linewidth=3)

						lines!(Fig[1,2], [Point(log1p(0.0), log1p(KsMat[iZ] * cst.MmS_2_MmH)), Point(log1p(cst.Mm_2_kPa * hydro.ΨmacMat[iZ]), log1p.(KsMat[iZ]* cst.MmS_2_MmH))], color=:brown, linewidth=3)
						

					# TAGGING
						if option.data.Kθ
							Leg = Fig[1, end+1] = Legend(Fig, [Fig_θΨobs, Fig_θΨsim, Fig_TotalPorosity, Fig_Kθobs, Fig_Kθsim, Fig_Ks], ["θobs(Ψ)", "θsim(Ψ)", "Φ", "Kobs(Ψ)", "Ksim(Ψ)", "Ksₛᵢₘ"])
						else
							Leg = Fig[1, end+1] = Legend(Fig, [Fig_θΨobs, Fig_θΨsim, Fig_TotalPorosity, Fig_Kθsim, Fig_Ks], ["θobs(Ψ)", "θsim(Ψ)", "Φ", "Ksim(Ψ)", "Ksₛᵢₘ"])
						end

					Fig[2, 1:2] = Leg
					trim!(Fig.layout)
					Leg.orientation = :horizontal
					Leg.tellheight = true
					
					Path = path.plotSoilwater.Plot_θΨK * "Lab_ThetaH_" * string(path.option.ModelName) * "_" * string(IdSelect[iZ]) * ".svg" 
					save(Path, Fig)
	
					# Displaying figure in VScode
					if option.general.PlotVscode
						display(Fig)
					end
				
				end # for iZ
				
			# ------------------------END: Plotting---------------------------  
			println("  ==  END: Plotting HydroParam  == \n")		
			return nothing
			end  # function: HYDROPARAM
	
	end  # module lab
	# ............................................................


	# =============================================================
	#		module: ksmodel

	# =============================================================
	module ksmodel
		import ...cst, ...θψ_2_KsψModel
		using CairoMakie, ColorSchemes
		import Polynomials
		import SpecialFunctions: erfc, erfcinv

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KSMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSMODEL(KₛModel, KΨ_Obs₁₀ₖₚₐ, KΨ_Sim₁₀ₖₚₐ, Ksₒᵦₛ, NameSim::String, Path::String, θrₒᵦₛ, θsₒᵦₛ, σₒᵦₛ, option)

			# Title
				Title = "K(Ψ)model" * option.ksModel.KₛModel⍰[end-1:end]

			# Dimension of figure
            Height = 1000 # Height of plot
            Width  = 1000  # Width of plot

			# Size of X and Y label
            XlabelSize = 45
            YlabelSize = 45
            NumberSize = 40

			# Title
				TitleSize = 60

			# Labels size of colourbar
             TickLabelSize = 35
             TickSize      = 20
             LabelSize     = 35
        
			# Colour map
				ColourMap = :plasma # :plasma, :ice, :viridis, :plasma

			# Activating the figure
				CairoMakie.activate!(type = "svg")
				Fig = Figure(font="Sans", fontsize=NumberSize)

			# PLOTTING KS	
			Axis_Ks = Axis(Fig[1,1], width=Width, height=Height, aspect=1, xlabel=L"$Ks _{sim}$ $[mm$ $h^{-1}]$", ylabel=L"$Ks _{model}$ $[mm$ $h^{-1}]$", xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false, xscale=Makie.pseudolog10, yscale=Makie.pseudolog10, xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9))

            Ksₒᵦₛ   = Ksₒᵦₛ .* cst.MmS_2_MmH
            KₛModel = KₛModel .* cst.MmS_2_MmH

            Ks_Min = minimum([minimum(Ksₒᵦₛ), minimum(KₛModel)])
            Ks_Max = maximum([maximum(Ksₒᵦₛ), maximum(KₛModel)])

				@show Ks_Min, Ks_Max

				# Ks_Max = 0.099371778 # mm/s
				
				xlims!(Axis_Ks, 0, Ks_Max)
				ylims!(Axis_Ks, 0, Ks_Max)

				# KsTicks = (range(0.0, stop=Ks_Max, length=10)) 
				Axis_Ks.xticks = [0, 1, 10, 10^2, 10^3, 10^4] 
				# (KsTicks, string.( floor.(KsTicks, digits=1)))
				# Axis_Ks.yticks = (KsTicks, string.(floor.(KsTicks, digits=1)))
				Axis_Ks.yticks =  [0, 1, 10, 10^2, 10^3, 10^4] 
				# Axis_Ks.xticklabelrotation = π/3

				ΔΘsMacΘr = θsₒᵦₛ .-  θrₒᵦₛ

				Fig_Ks = scatter!(Axis_Ks, Ksₒᵦₛ, KₛModel, color=σₒᵦₛ, markersize=135.0*ΔΘsMacΘr, marker=:circle, colormap=ColourMap, strokecolor=:black, strokewidth=1)
				Line = range(0.0, stop=Ks_Max, length=10) 
				Fig_Ks = lines!(Fig[1,1], Line, Line, color=:grey, linestyle=:dash, linewidth=5)

				# Leg1 = Colorbar(Fig, Fig_Ks, label = "Theta", ticklabelsize = 14, labelpadding = 5, width = 10)

			# PLOTTING K₁₀ₖₚₐ
				Axis_KΨ = Axis(Fig[1,2], aspect = 1, width= Width, height=Height, xlabel=L"$K10kpa _{sim}$ $[mm$ $h^{-1}]$", ylabel=L"$K10kpa _{model}$ $[mm$ $h^{-1}]$", xlabelsize=XlabelSize, ylabelsize=YlabelSize, xminorticksvisible=true, xminorticks=IntervalsBetween(9), yminorticksvisible=true, yminorticks=IntervalsBetween(9))

				KΨ_Obs₁₀ₖₚₐ = KΨ_Obs₁₀ₖₚₐ .* cst.MmS_2_MmH
				KΨ_Sim₁₀ₖₚₐ = KΨ_Sim₁₀ₖₚₐ .* cst.MmS_2_MmH

				KΨ_Obs₁₀ₖₚₐ_Min = minimum([minimum(KΨ_Sim₁₀ₖₚₐ), minimum(KΨ_Obs₁₀ₖₚₐ)])
				KΨ_Sim₁₀ₖₚₐ_Max = maximum([maximum(KΨ_Sim₁₀ₖₚₐ), maximum(KΨ_Obs₁₀ₖₚₐ)])

				KΨ_Sim₁₀ₖₚₐ_Max = 1.5

				xlims!(Axis_KΨ, 0.0, KΨ_Sim₁₀ₖₚₐ_Max)
				ylims!(Axis_KΨ, 0.0, KΨ_Sim₁₀ₖₚₐ_Max)

				Axis_KΨ.xticks = [0, 0.5, 1, 1.5, 2] 
				Axis_KΨ.yticks = [0, 0.5, 1, 1.5, 2] 

				# Axis_KΨ.xticklabelrotation = π/3

				Fig_KΨ = scatter!(Fig[1,2], KΨ_Obs₁₀ₖₚₐ, KΨ_Sim₁₀ₖₚₐ, color=σₒᵦₛ, markersize=135.0*ΔΘsMacΘr, marker=:circle, colormap=ColourMap, strokecolor=:black, strokewidth=1)

				Line = range(0.0, stop=Ks_Max, length=10) 
				Fig_Ks = lines!(Fig[1,2], Line, Line, color=:grey, linestyle=:dash, linewidth=5)
					
			# Colour bas
				Colorbar(Fig[1,3], limits=(minimum(σₒᵦₛ), maximum(σₒᵦₛ)+0.001), colormap =ColourMap, label="σ[-]", vertical=true, labelsize=LabelSize, width=30, ticksize=TickSize, ticklabelsize=TickLabelSize, labelpadding=5) # :thermal, :ice, :viridis, :plasma
				
			# Letters
				for (ax, label) in zip([Axis_Ks, Axis_KΨ], ["(A)", "(B)"])
					text!(
						ax, 0, 1,
						text = label, 
						font = :bold,
						align = (:left, :top),
						offset = (4, -2),
						space = :relative,
						fontsize = TitleSize,
						colour=:brown
					)
				end

			# Final adjustments
				Label(Fig[1, 1:2, Top()], Title, valign=:bottom, font=:bold, padding = (0, 0, 20, 0), fontsize=TitleSize)

				resize_to_layout!(Fig)
				trim!(Fig.layout)
				colgap!(Fig.layout, 20)
				rowgap!(Fig.layout, 20)

			Pathₛ = Path * "_" * NameSim * ".svg" 

			save(Pathₛ, Fig)
			# Displaying figure in VScode
				if option.general.PlotVscode
					display(Fig)
				end

		return nothing
		end  # function: KSMODEL
		# ------------------------------------------------------------------

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSMODEL_TCLAY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSMODEL_FUNCTIONS(Path, option, ksmodelτ, ipClass; τclayₘₐₓ=ksmodelτ.τclayₘₐₓ[ipClass], τclay₀=ksmodelτ.τclay₀[ipClass],τ₁ₐ=ksmodelτ.τ₁ₐ[ipClass], τ₂ₐ=ksmodelτ.τ₂ₐ[ipClass],τ₃ₐ=ksmodelτ.τ₃ₐ[ipClass], τclayΔθsr=ksmodelτ.τclayΔθsr[ipClass])
		
			# DERIVING THE DATA TO PLOT
				σ = 0.75:0.001:3.0
				Nσ  = length(σ)

				ΘsMacMatΘr =0.1:0.1:0.6
				NΘsΘr = length(ΘsMacMatΘr)
				Func_Ks1=fill(0.0, (NΘsΘr, Nσ))
				Func_Ks3=fill(0.0, (NΘsΘr, Nσ))
				ΨmacMat = 100.0
				T2_Max = 3.0; T3_Max = 4.0
				for iΘsΘr=1:NΘsΘr
					for iσ =1:Nσ
						T1 =  10.0 ^ (τ₁ₐ / (τ₁ₐ - 1.0))
						T2 = T2_Max * (1.0 - τ₂ₐ)
						T3 = T3_Max * (1.0 - τ₃ₐ)

						ΨmMean = exp((log(√ΨmacMat) + log(ΨmacMat)) * 0.5)
						Ψm = ΨmMean * exp(σ[iσ] * 3.0)

					# Ks model 1 not corrected for clay
						Func_Ks1[iΘsΘr, iσ] = 60.0 * 60.0 * T1 * cst.KunsatModel * π * ((ΘsMacMatΘr[iΘsΘr]) * ((cst.Y / Ψm) ^ T2) * exp(((T2 * σ[iσ]) ^ 2.0) / 2.0)) ^ T3

					# Ks model 2 corrected for clay
						Ψ_Clay =  160000.0 * ( ( (cst.Y  / 0.002) - (cst.Y / 0.5) ) / ((cst.Y  / 0.001) - (cst.Y  / 0.5)) ) ^ 2.0

						Clay = 0.5 * erfc((log(Ψ_Clay / Ψm)) / (σ[iσ] * √2.0))

						X_Clay₁ =  τclay₀

						Clayₙ = max(Clay - X_Clay₁, 0.0) / (1.0 - X_Clay₁)

						ΔθsMacθrₙ =  max(ΘsMacMatΘr[iΘsΘr] - τclayΔθsr , 0.0) 
						Tclay_Max =  1.0 + ΔθsMacθrₙ * (τclayₘₐₓ - 1.0) 

						Tclay = Tclay_Max - (Tclay_Max - 1.0) * cos(Clayₙ * π * 0.5) 

						Func_Ks3[iΘsΘr, iσ] = 60.0 * 60.0 * T1 * cst.KunsatModel * π * ((ΘsMacMatΘr[iΘsΘr]) ^ (Tclay / T3) * ((cst.Y / Ψm) ^ T2) * exp(((T2 * σ[iσ]) ^ 2.0) / 2.0)) ^ T3
					end
				end

			# PLOTTING
				#P arameters
				# Dimensions of figure
					Height = 800 # Height of plot
					Width  = 1000  # Width of plot

				# Size of X and Y label
					XlabelSize = 45
					YlabelSize = 45
					NumberSize = 40

				# Title size
					TitleSize = 60

				# Labels size of colourbar
					TickLabelSize = 35
					TickSize      = 20
					LabelSize     = 30
			
				# Colour map
					ColourMap = :plasma # :thermal :plasma, :ice, :viridis, :plasma

			# Activating the figure
				CairoMakie.activate!(type = "svg")
				Fig = Figure(font="Sans", fontsize=NumberSize)

			# PLOTTING KsModel1	
				Axis_KsModel1 = Axis(Fig[1,1], title="", width=Width, height=Height, xlabel=L"$σ [-]$", ylabel=L"$Ks _{model}$ $[mm$ $h^{-1}]$", xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false, xminorticksvisible=true, xminorticks=IntervalsBetween(10), yminorticksvisible=true, yminorticks=IntervalsBetween(10), yscale=Makie.pseudolog10)

				Axis_KsModel1.yticks =  [0, 10^0, 10^1, 10^2, 10^3, 10^4] 

				Colormap = cgrad(colorschemes[ColourMap], NΘsΘr, categorical = true)
				for iΘsΘr=1:NΘsΘr
					Fig_Model1 = lines!(Axis_KsModel1, σ[1:Nσ], Func_Ks1[iΘsΘr, 1:Nσ], linewidth=5, colormap =Colormap[iΘsΘr], label =string(ΘsMacMatΘr[iΘsΘr]))
				end

			# PLOTTING KsModel2	
				Axis_KsModel3 = Axis(Fig[2,1], title="", width=Width, height=Height, xlabel=L"$σ [-]$", ylabel=L"$Ks _{model}$ $[mm$ $h^{-1}]$", xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false, xminorticksvisible=true, xminorticks=IntervalsBetween(10), yminorticksvisible=true, yminorticks=IntervalsBetween(10), yscale=Makie.pseudolog10)

				Axis_KsModel3.yticks =  [0, 10^0, 10^1, 10^2, 10^3, 10^4] 

				Colormap = cgrad(colorschemes[:thermal], NΘsΘr, categorical = true)
				for iΘsΘr=1:NΘsΘr
					Fig_Model3 = lines!(Axis_KsModel3, σ[1:Nσ], Func_Ks3[iΘsΘr, 1:Nσ], linewidth=5, colormap =Colormap[iΘsΘr])
				end
				
				Leg = Legend(Fig[1:2,2], Axis_KsModel1, "θₛ-θᵣ", framevisible=true, tellheight=true, tellwidth=true, labelsize=LabelSize, margin=(30, 30, 30, 30))

				# Letters
					for (ax, Label) in zip([Axis_KsModel1, Axis_KsModel3], ["(i)", "(ii)"])
						text!(
							ax, 3.0, 150,
							text = Label, 
							font = :bold,
							align = (:right, :top),
							# offset = (40, -2),
							# space = :relative,
							fontsize = 40,
						)
					end

			# Final adjustments
				resize_to_layout!(Fig)
				trim!(Fig.layout)
				colgap!(Fig.layout, 20)
				rowgap!(Fig.layout, 20)

			Pathₛ = Path * "_" * "Func_KsModel" * ".svg" 

			save(Pathₛ, Fig)
			# Displaying figure in VScode
				if option.general.PlotVscode
					display(Fig)
				end
		return nothing
		end  # function: KSMODEL_TCLAY
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSMODEL_TCLAY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSMODEL_TCLAY(Path, option, ksmodelτ, ipClass; τclayₘₐₓ=ksmodelτ.τclayₘₐₓ[ipClass], τclay₀=ksmodelτ.τclay₀[ipClass], τclayΔθsr=ksmodelτ. τclayΔθsr[ipClass])
			# DERIVING THE DATA TO PLOT
				Tclay_Min = 1.0
				
				X_Clay₁ =  τclay₀ # τclay₀
				Clay = 0.0:0.001:1.0
				Nclay = length(Clay)

				ΘsΘr =0.1:0.1:0.6
				NΘsΘr = length(ΘsΘr)

				Func_Tclay=fill(0.0, (NΘsΘr, Nclay))
 				
				for iΘsΘr=1:NΘsΘr
					for iClay =1:Nclay

						X_Clay₁ =  τclay₀

						Clayₙ = max(Clay[iClay] - X_Clay₁, 0.0) / (1.0 - X_Clay₁)

						ΔθsMacθrₙ =  max(ΘsΘr[iΘsΘr] - τclayΔθsr , 0.0) 

						Tclay_Max =  1.0 + ΔθsMacθrₙ * (τclayₘₐₓ - 1.0) 

						Tclay = Tclay_Max - (Tclay_Max - 1.0) * cos(Clayₙ * π * 0.5) 

						Func_Tclay[iΘsΘr, iClay] =	ΘsΘr[iΘsΘr] ^ Tclay	
					end
				end

			# PLOTTING
				#P arameters
				# Dimensions of figure
					Height = 800 # Height of plot
					Width  = 1000  # Width of plot

				# Size of X and Y label
					XlabelSize = 45
					YlabelSize = 45
					NumberSize = 40

				# Title size
					TitleSize = 60

				# Labels size of colourbar
					TickLabelSize = 35
					TickSize      = 20
					LabelSize     = 35
			
				# Colour map
					ColourMap = :viridis # :plasma, :ice, :viridis, :plasma

			# Activating the figure
				CairoMakie.activate!(type = "svg")
				Fig = Figure(font="Sans", fontsize=NumberSize)

			# PLOTTING Tclay	
				Axis_Tclay = Axis(Fig[1,1], width=Width, height=Height, xlabel=L"$Clay$", ylabel=L"$ (θ_{s} - θ_{r}) ^{T_{clay}}$", xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false, xminorticksvisible=true, xminorticks=IntervalsBetween(10), yminorticksvisible=true, yminorticks=IntervalsBetween(10), yscale=log10, xlabelpadding=30)

				Axis_Tclay.xticks = [0, 0.25, 0.5, 0.75, 1] 
				xlims!(Axis_Tclay, 0, Clay[Nclay])
				ylims!(Axis_Tclay, 10^-6, ΘsΘr[NΘsΘr] )

				Fig_Tclay = empty
				for iΘsΘr=1:NΘsΘr
					Colormap = cgrad(colorschemes[ColourMap], NΘsΘr, categorical = true)
					Fig_Tclay = lines!(Axis_Tclay, Clay[1:Nclay], Func_Tclay[iΘsΘr, 1:Nclay], linewidth=6, colormap=Colormap[iΘsΘr], label =string(ΘsΘr[iΘsΘr]))
				end
				Leg = Legend(Fig[1,2], Axis_Tclay, "θₛ-θᵣ", framevisible=true, tellheight=true, tellwidth=true, labelsize=40, margin=(30, 30, 30, 30))
				# axislegend(Fig[1,2]; nbanks = 3, framecolor = (:grey, 0.5));

			# Final adjustments
				resize_to_layout!(Fig)
				trim!(Fig.layout)
				colgap!(Fig.layout, 20)
				rowgap!(Fig.layout, 20)

			Pathₛ = Path * "_" * "Tclay" * ".svg" 

			save(Pathₛ, Fig)
			# Displaying figure in VScode
				if option.general.PlotVscode
					display(Fig)
				end
		return nothing
		end  # function: KSMODEL_TCLAY
		# ------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KSMODEL_RF
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function KSMODEL_RF(Path, hydro, option, ksmodelτ, ipClass;τ₁ₐ=ksmodelτ.τ₁ₐ[ipClass], τ₂ₐ=ksmodelτ.τ₂ₐ[ipClass],τ₃ₐ=ksmodelτ.τ₃ₐ[ipClass],τ₂ₐMac=ksmodelτ.τ₂ₐMac[ipClass], τ₃ₐMac=ksmodelτ.τ₃ₐMac[ipClass])

			
			""" The rock corection is already performed in θ(Ψ) and therefore Ks is already corected. Nevertheles, the model is wrong for RF > Rf_StartIncrease as the Ks starts to increase again"""
					function ROCKCORRECTION!(hydro, iZ, Rf, θr, θs, θsMacMat; Rf_StartIncrease=0.4, Rf_EndIncrease=0.9, θs_Amplify=1.)

						Rf = min(Rf, Rf_EndIncrease)

						if Rf > Rf_StartIncrease
							# X values
								X = [Rf_StartIncrease, Rf_EndIncrease]
		
							# θs ----
								θs_NoRf = θs * 1.0
								Y_θs = [ (1.0 - Rf_StartIncrease) * θs_NoRf, θs_Amplify * θs_NoRf]
								Fit_θs = Polynomials.fit(X, Y_θs, 1)
								# θs = max(min(Fit_θs(Rf), hydro.θs_Max[iZ]), hydro.θs_Min[iZ])
								θs = Fit_θs(Rf)

							# θr ----
								θr_NoRf = θr * 1.0
								Y_θr = [(1.0 - Rf_StartIncrease) * θr_NoRf, θr_NoRf]
								Fit_θr = Polynomials.fit(X, Y_θr, 1)
								θr = max(min(Fit_θr(Rf), hydro.θr_Max[iZ]), hydro.θr_Min[iZ])

							# θsMacMat ----
								θsMacMat_NoRf =  θsMacMat * 1.0
								Y_θsMacMat = [min((1.0 - Rf_StartIncrease) * θsMacMat_NoRf, θs), 0.7 * (θs - θr) + θr]
								Fit_θsMacMat = Polynomials.fit(X, Y_θsMacMat, 1)	
								θsMacMat = min(Fit_θsMacMat(Rf), θs)

						else
							θs = θs * (1.0 - Rf)
							θr = θr * (1.0 - Rf)
							θsMacMat =  θsMacMat * (1.0 - Rf)
						end
					return θr, θs, θsMacMat
					end  # function: ROCKCORRECTION
				# ------------------------------------------------------------------

				# DERIVING THE DATA TO PLOT
					T2_Max = 3.0; T3_Max = 4.0
					T1     = 10.0 ^ (τ₁ₐ / (τ₁ₐ - 1.0))
					T2     = T2_Max * (1.0 - τ₂ₐ)
				
					T3     = T3_Max * (1.0 - τ₃ₐ)
					T1Mac  = T1
					T2Mac  = T2_Max * (1.0 - τ₂ₐMac)
					T3Mac  = T3_Max * (1.0 - τ₃ₐMac)

					ΨmacMat = 100.0
					Ψ₁ = 0.0
					σMac = hydro.σMac[1]
					ΨmMac = hydro.ΨmMac[1]
					ΨmMean = exp((log(√ΨmacMat) + log(ΨmacMat)) * 0.5)

					θr = [0.0,  0.0, 0.0, 0.0, 0.0]
					σ = [0.7, 0.8, 1.0, 1.5, 3.0]
					θs = [0.4, 0.4, 0.4, 0.4, 0.4]
					θsMacMat = θs .* 0.8
					Nsoil = length(θs)

					Rf = collect(0.0:0.001:0.9)
					Nrf = length(Rf)
							

					KsModel = fill(0.0::Float64, Nsoil, Nrf)
					for iSoil=1:Nsoil 
						for iRf=1:Nrf 							
							Ψm = ΨmMean * exp(σ[iSoil] * 3.0)

							θr₀, θs₀, θsMacMat₀ =  ROCKCORRECTION!(hydro, 1, Rf[iRf], θr[iSoil], θs[iSoil], θsMacMat[iSoil])

							KsModel[iSoil, iRf] =  60.0 * 60.0 * θψ_2_KsψModel.KsΨMODEL_NOINTEGRAL(T1, T1Mac, T2, T2Mac, T3, T3Mac, θr₀, θs₀, θsMacMat₀, σ[iSoil], σMac, Ψ₁, Ψm, ΨmMac)	
						end 
					end
					
				# PLOTTING
					# Dimensions of figure
						Height = 800 # Height of plot
						Width  = 1000  # Width of plot

					# Size of X and Y label
						XlabelSize = 45
						YlabelSize = 45
						NumberSize = 40

					# Title size
						TitleSize = 60

					# Labels size of colourbar
						TickLabelSize = 35
						TickSize      = 20
						LabelSize     = 30
				
					# Colour map
						ColourMap = :plasma # :thermal :plasma, :ice, :viridis, :plasma

					# Activating the figure
						CairoMakie.activate!(type = "svg")
						Fig = Figure(font="Sans", fontsize=NumberSize)

				# PLOTTING KsModel1	
					Axis_KsModel1 = Axis(Fig[1,1], title="", width=Width, height=Height, xlabel=L"$Rock Fragments$ $[%]$", ylabel=L"$Ks _{model}$ $[mm$ $h^{-1}]$", xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false,  yminorticksvisible=true, yminorticks=IntervalsBetween(10))

					Axis_KsModel1.xticks =  [0, 0.2, 0.4, 0.6, 0.8, 1] 
					Axis_KsModel1.yticks =  [0, 10, 20, 30, 40, 50, 60] 

					Colormap = cgrad(colorschemes[ColourMap], Nsoil, categorical = true)
					for iSoil=1:Nsoil
						Fig_Model1 = lines!(Axis_KsModel1, Rf, KsModel[iSoil, :], linewidth=5, colormap =Colormap[iSoil], label =string(σ[iSoil]))
					end
		
				Leg = Legend(Fig[1:2,2], Axis_KsModel1, "σ", framevisible=true, tellheight=true, tellwidth=true, labelsize=LabelSize, margin=(30, 30, 30, 30))


			# Final adjustments
				resize_to_layout!(Fig)
				trim!(Fig.layout)
				colgap!(Fig.layout, 20)
				rowgap!(Fig.layout, 20)

			Pathₛ = Path * "_" * "Func_RockFragment" * ".svg" 

			save(Pathₛ, Fig)
			# Displaying figure in VScode
				if option.general.PlotVscode
					display(Fig)
				end
			return nothing	
			end  # function: KSMODEL_RF
		# ------------------------------------------------------------------


	end  # module: ksmodel

	
	# ............................................................

	# =============================================================
	#		MODULE: psd
	# =============================================================
	module psd
		using Plots, Plots.PlotMeasures, LaTeXStrings
		import ...wrc, ...kunsat, ...cst, ...psdThetar, ...psdFunc, ...bestFunc
		export PLOT_θr, PLOT_IMP_MODEL, PLOT_PSD_θΨ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_θr
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PLOT_θr(∑Psd, hydro, hydroPsd, NiZ, param, Path)
				println("  ==  START: Plotting PLOT_θr  ==")

				# Sorting ascending order with clay fraction
					Array      = zeros(Float64, 3, length(∑Psd[1:NiZ, param.psd.Psd_2_θr_Size]))
					Array      = zeros(Float64, (3, NiZ))
				
					Array[1,:] = ∑Psd[1:NiZ, param.psd.Psd_2_θr_Size] # Clay fraction
					Array[2,:] = hydroPsd.θr[1:NiZ]
					Array[3,:] = hydro.θr[1:NiZ]
					Array      = sortslices(Array, dims=2)
					Clay       = Array[1,:] # Clay fraction
					θr_Psd     = Array[2,:]
					θr         = Array[3,:]
				
				# Minimum and maximum value
					θr_Min = 0.01 

					θr_Max = maximum(hydroPsd.θr_Max) + 0.05
					Clay_Min = 0.1
					Clay_Max = maximum(∑Psd[1:NiZ, param.psd.Psd_2_θr_Size]) + 0.05
				
				# PLOT 1 <>=<>=<>=<>=<>=<>
					# pgfplotsx()
					# Plot θr(Clay)
						X = Clay
						Y = θr
						Plot_θr = Plots.plot(X, Y, seriestype=:scatter, label=L"\theta _{r}", color= :violet, shape= :square, markersize=4, legend=:bottomright, size=(5000,400))

					# Plot θr_psd(Clay)
						X = Clay
						Y = θr_Psd
						Plots.plot!(X ,Y, seriestype=:line, label=L"\theta _{r psd}", color= :blue, lw=2)
				
					# General attributes
						xlabel!(L"Clay \ [g \ g^{-1}]")                         
						ylabel!(L"\theta _{r} [cm^{3} \ cm^{-3}]")
						Plots.plot!(xlims= (Clay_Min, Clay_Max), ylims= (θr_Min, θr_Max))

				# PLOT 2 <>=<>=<>=<>=<>=<>
					# Plot θr_Psd(θr)
						X = θr
						Y = θr_Psd
						Plot_θr_Psd = Plots.plot(X ,Y, seriestype=:scatter, color=:violet, shape=:square, markersize=4, size=(800,400))
						
					# 1:1 line
						X = range(θr_Min, stop=θr_Max, length=10)
						Y = X
						Label = "1:1"
						Plots.plot!(X, Y, seriestype=:line, label= Label , color= :black, linestyle= :dot, lw=2)

					# General attributes
						xlabel!(L"\theta _{r} [cm^3 cm^{-3}]")
						ylabel!(L"\theta _{r \ psd} [cm^3 cm^{-3}]")
						Plots.plot!(xlims= (θr_Min, θr_Max), ylims= (θr_Min, θr_Max))

				Plot = Plots.plot(Plot_θr, Plot_θr_Psd)
				Plots.savefig(Plot, Path)
				println("    ~  $(Path) ~")
			
			println("  ==  END: Plotting PLOT_θr  == \n")
			return nothing
			end # function: PLOT_θr


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_IMP_MODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PLOT_IMP_MODEL(∑Psd, hydro, IdSelect, NiZ, N_Psd, option, param, Path, Psd, Rpart)
				println("  ==  START: PLOT_IMP_MODEL  ==")	

				for iZ = param.globalparam.N_iZ_Plot_Start: param.globalparam.N_iZ_Plot_End
					Rpart_Min = minimum(Rpart[iZ,1:N_Psd[iZ]])
					Rpart_Max = maximum(Rpart[iZ,1:N_Psd[iZ]])

					∑Psd_Min  = minimum(∑Psd[iZ,1:N_Psd[iZ]])
					∑Psd_Max  = maximum(∑Psd[iZ,1:N_Psd[iZ]])

					Psd_Min  = minimum(Psd[iZ,1:N_Psd[iZ]])
					Psd_Max  = maximum(Psd[iZ,1:N_Psd[iZ]])

					IntergranularMixing = zeros(Float64, N_Psd[iZ])
					ξ = zeros(Float64, N_Psd[iZ])
					for iRpart = 1:N_Psd[iZ]
						# added at a later stage
						ξ2 = psdFunc.imp.∑PSD_2_ξ2(∑Psd[param.psd.imp.∑Psd_2_ξ2_Size], param)

						ξ[iRpart] = psdFunc.imp.INTERGRANULARMIXING(param, Rpart[iZ,iRpart], param.psd.imp.ξ1, ξ2)

						IntergranularMixing[iRpart] = (Rpart[iZ, iRpart] ^ -ξ[iRpart]) 
					end # for iRpart = 1:N_Psd[iZ]

					# << PLOT 1 >>
						# Plot_∑Psd_Rpart
							X = Rpart[iZ,1:N_Psd[iZ]]
							Y = ∑Psd[iZ,1:N_Psd[iZ]]
							Plot_∑Psd_Rpart = Plots.plot(X ,Y, seriestype=:scatter, color= :teal, shape= :square, markersize= 4, size=(800,400))
							Plots.plot!(X ,Y, seriestype=:line, color= :teal)

						# Plot_∑Psd_Rpart: General attributes
							xlabel!(L"R_{part} [mm]")
							ylabel!(L"\sum \ PSD")
							Plots.plot!(xlims= (Rpart_Min, Rpart_Max), ylims= (∑Psd_Min, ∑Psd_Max), xscale= :log10)

					# << PLOT 2 >>
						# Plot_Psd_Rpart
							X = Rpart[iZ,1:N_Psd[iZ]]
							Y = Psd[iZ,1:N_Psd[iZ]]
							# Plot_Rpart_Psd = Plots.plot(X ,Y, seriestype=:scatter, color= :blue, shape= :square, markersize= 4, size=(800,400))
							Plot_Rpart_Psd = Plots.plot(X ,Y, seriestype=:line, color= :blue, shape= :square, markersize= 4, size=(800,400))

							xlabel!(L"R_{part} \ [mm]")
							ylabel!(L"PSD [mm]")
							Plots.plot!(xlims= (Rpart_Min, Rpart_Max), ylims= (Psd_Min, Psd_Max), xscale= :log10)


					# << PLOT 3 >>
						# Plot NormMixing_Rpart
							X = Rpart[iZ,1:N_Psd[iZ]]
							Y = IntergranularMixing[1:N_Psd[iZ]] / maximum( IntergranularMixing[1:N_Psd[iZ]] )
							Plot_NormMixing_Rpart = Plots.plot(X, Y, seriestype=:line, color= :green)

							xlabel!(L"R_{part} \ [mm]")
							ylabel!(L"R_{part}^{-\xi(R_{Part})}")
							Plots.plot!(xlims= (Rpart_Min, Rpart_Max), ylims= (0.0, 1.1), xscale= :log10)

					Plot = Plots.plot(Plot_∑Psd_Rpart, Plot_Rpart_Psd, Plot_NormMixing_Rpart, layout = (3,1))
					Path₀ = Path * "IMP_" * string(option.hydro.HydroModel⍰) * "_" *string(IdSelect[iZ]) * ".svg"
					Plots.savefig(Plot, Path₀)
					# println("    ~  $(Path₀) ~")
				end # for iZ
			println("  ==  END: PLOT_IMP_MODEL  == \n")
			return nothing	
			end # function: PLOT_IMP_MODEL


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_IMP_ΘΨ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PLOT_PSD_θΨ(hydro, hydroPsd, IdSelect, NiZ, N_Psd, N_θΨobs, option, param, Path, θ_Rpart, θ_θΨobs, Ψ_Rpart, Ψ_θΨobs; N_Se= 100)
			
				println("  ==  START: Plotting PLOT_PSD_θΨ  ==")

				θ_θΨobs_Psd = fill(0.0::Float64, (N_Se))

				for iZ = param.globalparam.N_iZ_Plot_Start:param.globalparam.N_iZ_Plot_End
					# Range of plots
						Ψ_θΨobs_Min = 10.0 ^ 0.0 # [mm]

						Ψ_θΨobs_Max = 150000 + 100000 # [mm]

						Ψ_Sim = 10.0 .^ range(log(Ψ_θΨobs_Min), stop=log(Ψ_θΨobs_Max), length=N_Se)

						θ_θΨobs_Max = hydroPsd.Φ[iZ] + 0.1

					# Simulated 
						for iΨ = 1:N_Se
							θ_θΨobs_Psd[iΨ] = wrc.Ψ_2_θ(option.psd,Ψ_Sim[iΨ], iZ, hydroPsd)
						end # iΨ 

					# Plot_θ_Ψ: Psd model fitting for e.g. Kosugi model
						X = Ψ_Sim[1:N_Se] .* cst.Mm_2_Cm
						Y = θ_θΨobs_Psd[1:N_Se]
						Label = "PsdKosugi"
						Plot_θ_Ψ_Psd = Plots.plot(X ,Y, seriestype=:line, label=Label, color= :blue, lw=2)

					# Plot_θ_Ψ: Psd model points
						X = Ψ_Rpart[iZ, 1:N_Psd[iZ]] .* cst.Mm_2_Cm
						Y = θ_Rpart[iZ, 1:N_Psd[iZ]]
						Label = "PsdModel"
						Plot_θ_Ψ_Psd = Plots.plot!(X ,Y, seriestype=:scatter, label=Label, color= :violet, shape= :circle, markersize=4)

					# Plot_θ_Ψ: Observed
					if option.run.HydroLabθΨ⍰ ≠ "No" 
						X = Ψ_θΨobs[iZ,1:N_θΨobs[iZ]] .* cst.Mm_2_Cm
						Y = θ_θΨobs[iZ,1:N_θΨobs[iZ]]
						Label = "LabObs"
						Plot_θ_Ψ = Plots.plot!(X ,Y, seriestype=:scatter, label=Label, color= :red, shape= :square, markersize=4)

						# Plot_θ_Ψ: Total porosity point
							X = zeros(Float64,1)
							X[1] = Ψ_θΨobs_Min * cst.Mm_2_Cm
							Y = zeros(Float64,1)
							Y[1] = hydro.Φ[iZ]
							Label = "\$ \\phi \$"
							Plots.plot!(X, Y, seriestype=:scatter, label= Label, color= :green, shape= :square, markersize=4) 
					end

					# Plot_θ_Ψ: General attributes
						xlabel!(L"\psi \ [cm]")
						ylabel!(L"\theta \ [cm^3 cm^{-3}]")
						Plots.plot!(xlims =(Ψ_θΨobs_Min*cst.Mm_2_Cm, Ψ_θΨobs_Max*cst.Mm_2_Cm), ylims =(0.0, θ_θΨobs_Max), xscale= :log10, size=(800,400))

					Path₀ = Path * "Psd_ThetaH_" * string(option.hydro.HydroModel⍰) * "_" *string(IdSelect[iZ]) * ".svg"     
					Plot = Plots.plot(Plot_θ_Ψ_Psd)
					Plots.savefig(Plot, Path₀)
					# println("    ~  $(Path₀) ~")
				end # iZ
			println("  ==  END: Plotting PLOT_PSD_θΨ  == \n")
			return	nothing	
			end # function PLOT_IMP_ΘΨ

	end  # module: psd
	# ............................................................


	# =============================================================
	#		MODULE: infilt
	# =============================================================
	module infilt
		import ...wrc, ...kunsat, ...cst, ...psdThetar, ...psdFunc, ...bestFunc, ...sorptivity
		using Plots, Plots.PlotMeasures, LaTeXStrings
		export  PLOT_∑INFILT, PLOT_∑INFILT_θΨ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_∑INFILT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PLOT_∑INFILT(∑Infilt_1D, ∑Infilt_3D, ∑Infilt_Obs, IdSelect, infiltOutput, N_Infilt, NiZ, option, param, Path, Tinfilt)
			println("  ==  START: PLOT_∑INFILT  == \n")
			
				for iZ = param.globalparam.N_iZ_Plot_Start:param.globalparam.N_iZ_Plot_End
					# << PLOT 1 >>
						Title = " iZ= $(IdSelect[iZ])"
						# Plot_∑infilt_Obs

							Label ="Obs_$(string(option.infilt.DataSingleDoubleRing⍰))_Ring"
							X = Tinfilt[iZ,1:N_Infilt[iZ]] / 60.0
							Y = ∑Infilt_Obs[iZ,1:N_Infilt[iZ]]
							Plot_∑infilt_Obs = Plots.plot(X, Y, seriestype=:scatter, label=Label, color= :red, shape= :square, markersize=4, marker = (Plots.stroke(1, :red))) 

						# Plot_∑infilt_Sim
							Label = "Sim_3D"
							X = Tinfilt[iZ,1:N_Infilt[iZ]] / 60.0
							Y = ∑Infilt_3D[iZ,1:N_Infilt[iZ]]
							Plots.plot!(X, Y, seriestype=:line, label=Label, color= :blue, shape= :square, markersize=4, marker = (Plots.stroke(1, :blue))) 


							Label = "Sim_1D"
							Y2 = ∑Infilt_1D[iZ,1:N_Infilt[iZ]]
							Plots.plot!(X, Y2, seriestype=:line, label=Label, color= :green, shape= :square, markersize=4,  marker = (Plots.stroke(1, :green))) 

						# TransSteady
							Label="T_TransSteady"
							X = zeros(Float64,1)
							Y = zeros(Float64,1)
							X[1] = Tinfilt[iZ,infiltOutput.iT_TransSteady_Data[iZ]] / 60.0
							Y[1] = ∑Infilt_Obs[iZ,infiltOutput.iT_TransSteady_Data[iZ]]
							Plots.plot!(X, Y, seriestype=:scatter, label=Label, color= :violet, shape= :circle, markersize=10, title=Title) 

							Plots.xlabel!(L"Time [minutes]")
							Plots.ylabel!(L"\sum infiltration \ [mm]")      
							
						Path₂ = Path * "INFIL_" * string(option.infilt.Model⍰)  *  "_" * string(IdSelect[iZ]) *  ".svg"

					Plots.savefig(Plot_∑infilt_Obs, Path₂)
					# println("    ~  $(Path₂) ~")
				end # for iZ=1:NiZ
			println("  ==  END: PLOT_∑INFILT  == \n")
			return nothing
			end # PLOT_∑INFILT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PLOT_∑INFILT_θΨ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PLOT_∑INFILT_θΨ(hydroInfilt, IdSelect, NiZ, optim, option, param, Path; hydro=[], N_Se=100)
			println("  ==  START: PLOT_∑INFILT_θΨ  ==")

				θ_Infilt      = fill(0.0::Float64, (N_Se))
				θ_Obs         = fill(0.0::Float64, (N_Se))
				Kunsat_Infilt = fill(0.0::Float64, (N_Se))
				Kunsat_Obs    = fill(0.0::Float64, (N_Se))

				for iZ = param.globalparam.N_iZ_Plot_Start: param.globalparam.N_iZ_Plot_End	
					Ψ_θΨobs_Min = 10.0 ^ -2 # [mm]

					Ψ_θΨobs_Max = 200000.0 * 10.0 # [mm]

					Ψ = 10.0 .^ range(log(Ψ_θΨobs_Min), stop=log(Ψ_θΨobs_Max), length=N_Se)

					θ_θΨobs_Max = hydroInfilt.Φ[iZ] + 0.1

					if option.run.HydroLabθΨ⍰ ≠ "No" && "Ks" ∈ optim.ParamOpt
						K_Ψ_Max = max(hydroInfilt.Ks[iZ], hydro.Ks[iZ]) * 1.1
					else
						K_Ψ_Max = hydroInfilt.Ks[iZ] * 1.1
					end #  "Ks" ∈ optim.ParamOpt

					for iΨ = 1:N_Se
						θ_Infilt[iΨ] = wrc.Ψ_2_θ(option.infilt,Ψ[iΨ], iZ, hydroInfilt)

						Kunsat_Infilt[iΨ] = kunsat.KUNSAT_θΨSe(option.infilt, Ψ[iΨ], iZ, hydroInfilt)

						if option.run.HydroLabθΨ⍰ ≠ "No"
							θ_Obs[iΨ] = wrc.Ψ_2_θ(option.infilt,Ψ[iΨ], iZ, hydro)

							if option.run.HydroLabθΨ⍰ ≠ "No" && "Ks" ∈ optim.ParamOpt
								Kunsat_Obs[iΨ] = kunsat.KUNSAT_θΨSe(option.infilt, Ψ[iΨ], iZ, hydro)
							end # "Ks" ∈ optim.ParamOpt		
						end # option.run.HydroLabθΨ⍰ ≠ :No
					end # iΨ 

					#PLOT 1:  Plot_θ_Ψ
						# Plot_θ_Ψ: Simulated Infiltration
							X = Ψ[1:N_Se] .* cst.Mm_2_Cm
							Y = θ_Infilt[1:N_Se]
							Label = "Infiltration"
							Plot_θ_Ψ = Plots.plot(X, Y, seriestype=:line, label=Label, color= :blue, lw=2)

						# Plot_θ_Ψ: Observed
						if option.run.HydroLabθΨ⍰ ≠ "No"
							X = Ψ[1:N_Se] .* cst.Mm_2_Cm
							Y = θ_Obs[1:N_Se]
							Label = "Obs"
							Plot_θ_Ψ = Plots.plot!(X ,Y, seriestype=:line, label=Label, color= :red, lw=2)
						end # option.run.HydroLabθΨ⍰ ≠ :No

						# Plot_θ_Ψ: General attributes
							Plots.xlabel!("\\psi [cm]")
							Plots.ylabel!(L"\theta \ [cm^3 cm^{-3}]")
							Plots.plot!(xlims =(10.0*Ψ_θΨobs_Min*cst.Mm_2_Cm, Ψ_θΨobs_Max*cst.Mm_2_Cm), ylims =(0.0, θ_θΨobs_Max), xscale= :log10, size=(800,400), legend=:bottomleft)

						# PLOT2: Kunsat
							# Plot_K_Ψ: Obs K_Ψ
							X = Ψ[1:N_Se] .* cst.Mm_2_Cm
							Y = Kunsat_Infilt[1:N_Se] .* cst.MmS_2_CmH
							Label = "Infiltration"
							Plot_K_Ψ = Plots.plot(X, Y, seriestype=:line, label=Label, color= :blue, lw=2)

							# Plot_K_Ψ: Sim K_Ψ
							if option.run.HydroLabθΨ⍰ ≠ "No" && "Ks" ∈ optim.ParamOpt
								X = Ψ[1:N_Se] .* cst.Mm_2_Cm
								Y = Kunsat_Obs[1:N_Se] .* cst.MmS_2_CmH
								Label = "Obs"
								Plot_K_Ψ = Plots.plot!(X, Y, seriestype=:line, label=Label, color= :red, lw=2)
							end # "Ks" ∈ optim.ParamOpt

							# General attributes
								Plots.xlabel!("\\psi [cm]")
								Plots.ylabel!(L" K (\psi) \ [cm \ h^{-1}]")
								Plot_K_Ψ = Plots.plot!(xlims = (Ψ_θΨobs_Min*cst.Mm_2_Cm, Ψ_θΨobs_Max*cst.Mm_2_Cm), ylims = (10^-2.0, K_Ψ_Max * cst.MmS_2_CmH), xscale= :log10,  yscale= :log10, legend=:bottomleft, size=(800,400))

						Path₂ = Path * "Infilt_ThetaH_" * string(option.hydro.HydroModel⍰) * "_" *string(IdSelect[iZ]) * ".svg"     
						Plot = Plots.plot(Plot_θ_Ψ, Plot_K_Ψ)
						Plots.savefig(Plot, Path₂)
						# println("    ~  $(Path₂) ~")
				end # iZ

			println("  ==  END: PLOT_∑INFILT_θΨ  == \n")
			return nothing
			end  # function: PLOT_∑INFILT_θΨ

	end # module infilt
	# ............................................................

end  # module plot

