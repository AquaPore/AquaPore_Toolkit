# =============================================================
#		MODULE: hypixStart
# =============================================================

include("IncludeHypix.jl")

module hypixStart

	import ..cst, ..horizonLayer, ..hydroSmooth, ..hydroStruct, ..hypixModel, ..hypixOpt, ..memory, ..plotHypix, ..readHypix, ..readLinkingFile, ..stats, ..tableHypix, ..tool, ..waterBalance, ..Δtchange, ..θaver
	import Statistics: mean
	import Dates: now, value
	import CSV, Tables

	export HYPIX_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIX_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYPIX_START(NameLinkingFile, ProjectHypix)

		println("\n =============== Start running HyPix 1D =========================== \n")

		Time_Start = now()
		
		# GETTING PATHS ===
			Path_Hypix = dirname(dirname(@__DIR__)) # moving down the path twice

			dateHypix, Id, N_Scenario, pathInputHypix, SiteName = readLinkingFile.LINKING_FILE(Path_Hypix, NameLinkingFile, ProjectHypix)

		# READING VALUES FOR EVERY SCENARIOS ===
			∑∑ΔSink=[]; ∑ΔQ_Bot=[]; CccBest=[];  Efficiency=[];  Global_WaterBalance=[];  Global_WaterBalance_NormPr=[];  iNonConverge_iOpt=[];  NseBest=[];  θroot_Mean=[];  WilmotBest=[];  WofBest=[];  ΔRunTimeHypix=[];  ΔT_Average=[]; ∑Q_Z=[]; ∑Pr_Clim=[]; ∑Pet_Net=[]
		# SCENARIOS
			for iScenario = 1:N_Scenario #----------------------

				println("=== === === SITENAME=  ", ProjectHypix,"_", SiteName[iScenario], "  === === ===  \n")

				∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑PrThroughfall, ∑PrThroughfall_Climate, ∑T, ∑T_Climate, ∑T_Qobs, ∑T_Qobs, ∑ΔQ_Obs, ∑ΔQ_Obs, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, Flag_θΨini, Hpond, hydro_best, hydroHorizon, hydroHorizon_best, K_Aver_Vect, K_Aver₀_Vect, Layer, N_∑T_Climate, N_iRoot, N_Layer, N_SoilLayer, Nz, obsθ, optionHypix, paramHypix, pathInputHypix, pathOutputHypix, Pkₐᵥₑᵣ, Q, Residual, Temp, veg, veg_best, Z, ΔEvaporation, ΔLnΨmax, ΔPet, ΔPrThroughfall, ΔRootDensity, ΔRunoff, ΔSink, ΔT, θ, θini_or_Ψini, θSim, Ψ, Ψ_Min, Ψbest = readHypix.READ_START(dateHypix, Id, iScenario, N_Scenario, Path_Hypix, pathInputHypix, ProjectHypix, SiteName)

				# MEMORY FOR OUTPUT OF EVERY SCENARIO
					if iScenario == 1
						∑∑ΔSink, ∑Pet_Net, ∑Pr_Clim, ∑Q_Z, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iOpt, NseBest, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average, θroot_Mean = memory.MEMORY_SCENARIO(N_Scenario, optionHypix, paramHypix)
					end

				# OPTIMISATION
				# No optimisation
				if !(optionHypix.opt.Optimisation)
					paramHypix.opt.iOptMultiStep_Start = 1
					paramHypix.opt.iOptMultiStep_End   = 1
				end

				for iMultistep = paramHypix.opt.iOptMultiStep_Start:paramHypix.opt.iOptMultiStep_End

					# COUNT SIMULATIONS
						if optionHypix.opt.Optimisation

							iOpt_Count = iMultistep - paramHypix.opt.iOptMultiStep_Start + 1

							if iScenario ≥ 2
								error("HyPix error: if optimisation iScenario should be = 1")

							end

							println("		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
							println("		=== === === START: Multistep= ", iMultistep, " steps, \n")
						else
							iOpt_Count = iScenario
						end

					# AVERAGE INTERCELL_CONDUCTIVITY
						Pkₐᵥₑᵣ_Array=[1.0, 1.0, 1.0]
						Pkₐᵥₑᵣ[1] = Pkₐᵥₑᵣ_Array[1]
						for iZ=2:Nz-1
							Pkₐᵥₑᵣ[iZ] = Pkₐᵥₑᵣ_Array[2]
						end
						Pkₐᵥₑᵣ[Nz] = Pkₐᵥₑᵣ_Array[3]

					# OBTAINING HYDRAULIC AND VEGETATION PARAMETERS (depending of we have multistep optimisation)
					if optionHypix.opt.Optimisation

						hydroHorizon, optim, veg = readHypix.HYPIX_PARAM_OPT(hydroHorizon, iMultistep, N_SoilLayer, optionHypix, paramHypix, pathInputHypix.MultistepOpt[iScenario], veg)

						# Injecting hydraulic & vegetation parameters from output
						if optionHypix.opt.HydroVegParamReadFromOutput && paramHypix.opt.iOptMultiStep_Start == iMultistep
							@warn( "\n ========= HydroVegParamReadFromOutput = true ========= \n" )

							# Reading hydraulic parameters
								Path = pathOutputHypix.Table_Hydro  * "_" * string(paramHypix.opt.iOptMultiStep_Start-1) * ".csv"
								hydroHorizon, ~ = tool.readWrite.READ_STRUCT_SIMPLE(hydroHorizon, Path)
								
							# Reading vegetation parameters
								Path = pathOutputHypix.Table_Veg * "_" * string(paramHypix.opt.iOptMultiStep_Start-1) * ".csv"
								veg, ~ = tool.readWrite.READ_STRUCT_SIMPLE(veg, Path)
						end # INJECTING HYDRAULIC & VEGETATION PARAMETERS FROM OUTPUT
					else
						# options of optim		
							optim = (NparamOpt=0, Flag_Opt=false)		
					end # optionHypix.Optimisation

					@simd for iZ=1:N_SoilLayer
						hydroHorizon.So[iZ] = paramHypix.So # 1.0E-8
					end

					if !(optionHypix.opt.Optimisation) || optim.Flag_Opt == false
						if optionHypix.HydroSmooth
							hydro = hydroSmooth.HYDROHORIZON_2_HYDRO_SMOOTENING(discret, hydroHorizon, Layer, optionHypix)
						else
							hydro = horizonLayer.HYDROHORIZON_2_HYDRO(hydroHorizon, Layer, Nz, optionHypix)
						end
					else
						hydro = []
					end

					if optim.Flag_Opt
						hydro, hydro_best, hydroHorizon, hydroHorizon_best, veg, veg_best, WofBest = hypixOpt.HYPIXOPTIMISATION_START(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑PrThroughfall, ∑PrThroughfall_Climate, ∑T, ∑T_Climate, ∑T_Qobs, ∑ΔQ_Obs, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Flag_θΨini, Hpond, hydro_best, hydroHorizon, hydroHorizon_best, iOpt_Count, iScenario, K_Aver_Vect, K_Aver₀_Vect, Layer, N_∑T_Climate, N_iRoot, N_SoilLayer, Nz, obsθ, optim, optionHypix, paramHypix, pathInputHypix, Pkₐᵥₑᵣ, Q, Residual, veg, veg_best, WofBest, Z, ΔEvaporation, ΔLnΨmax, ΔPet, ΔPrThroughfall, ΔRootDensity, ΔRunoff, ΔSink, ΔT, θ, θini_or_Ψini, θSim, Ψ, Ψ_Min, Ψbest)
					end
				
					# if Flag_Opt then it will rerun with the optimal parameters
					∑Pet, ∑PrThroughfall, ∑T, ∑T_Climate, clim, discret, Hpond, iNonConverge, IterCount, N_iRoot, Nit, Nz, Q, veg, ΔEvaporation, ΔRootDensity, ΔRunoff, ΔT, θ, Ψ = hypixModel.HYPIX_MODEL(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet_Climate, ∑Pet, ∑PrThroughfall_Climate, ∑PrThroughfall, ∑T_Climate, ∑T, clim, CropCoeficientᵀ_η, CropCoeficientᵀ, discret, Flag_θΨini, Hpond, hydro, iScenario, K_Aver_Vect, K_Aver₀_Vect, N_∑T_Climate, N_iRoot, Nz, optionHypix, paramHypix, pathInputHypix, Pkₐᵥₑᵣ, Q, Residual, veg, Z, ΔEvaporation, ΔLnΨmax, ΔPet, ΔPrThroughfall, ΔRunoff, ΔSink, ΔRootDensity, ΔT, θ, θini_or_Ψini, Ψ_Min, Ψ, Ψbest)

					# WATER BALANCE
					# Computed after the warmup period
						∑∑WaterBalance, ∑WaterBalance_η, ∑ΔSink, i∑T_CalibrStart, ΔStorage = waterBalance.WATERBALANCE(∑T, clim, dateHypix, discret, hydro, N_iRoot, Nit, Nz, Q, ΔSink, ΔT, θ, Ψ)

						if optionHypix.θobs
							∑Runoff_Reduced, ∑T_Reduced, ∑WaterBalanceη_Reduced, ∑ΔQ_Obs_Reduced, Date_Reduced, Nit_Reduced, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPrGross_Reduced, ΔPrThroughfall_Reduced, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔSink_Reduced, ΔT_Reduced, θ_Reduced, θobs_Reduced, Ψ_Reduced = Δtchange.CHANGE_OUTPUT_ΔT(∑Pet[1:Nit], ∑PrThroughfall[1:Nit], ∑T[1:Nit], ∑T_Climate, ∑WaterBalance_η[1:Nit], ∑ΔSink[1:Nit], clim, dateHypix, Hpond[1:Nit], iScenario, Nit, Nz, optionHypix, paramHypix, pathInputHypix, Q[1:Nit,1:Nz+1], ΔEvaporation[1:Nit], ΔRunoff[1:Nit], ΔT[1:Nit], θ[1:Nit,1:Nz], Ψ[1:Nit,1:Nz], ∑ΔQ_Obs, ∑T_Qobs; obsθ=obsθ)							

						else
							∑Runoff_Reduced, ∑T_Reduced, ∑WaterBalanceη_Reduced, ∑ΔQ_Obs_Reduced, Date_Reduced, Nit_Reduced, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPrGross_Reduced, ΔPrThroughfall_Reduced, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔSink_Reduced, ΔT_Reduced, θ_Reduced, θobs_Reduced, Ψ_Reduced = Δtchange.CHANGE_OUTPUT_ΔT(∑Pet[1:Nit], ∑PrThroughfall[1:Nit], ∑T[1:Nit], ∑T_Climate, ∑WaterBalance_η[1:Nit], ∑ΔSink[1:Nit], clim, dateHypix, Hpond[1:Nit], iScenario, Nit, Nz, optionHypix, paramHypix, pathInputHypix, Q[1:Nit,1:Nz+1], ΔEvaporation[1:Nit], ΔRunoff[1:Nit], ΔT[1:Nit], θ[1:Nit,1:Nz], Ψ[1:Nit,1:Nz], ∑ΔQ_Obs, ∑T_Qobs)
						end

					# SUMMARY HOW GOOD THE SIMULATION
						# Computed climate day after the warmup period
							i∑T_CalibrStart_Day = 1::Int64 
							while ∑T_Climate[i∑T_CalibrStart_Day] < dateHypix.Δ∑T_StartSim
							# while ∑T_Climate[i∑T_CalibrStart_Day] < obsθ.∑T[1]
								i∑T_CalibrStart_Day += 1
							end
							i∑T_CalibrStart_Day += 1

						# Climate
							∑Pr_Clim[iOpt_Count] = sum(clim.Pr[i∑T_CalibrStart_Day:clim.N_Climate]) 
							∑Pet_Net[iOpt_Count] = sum(clim.Pet[i∑T_CalibrStart_Day:clim.N_Climate])

						# Soil water content for the rootzone at the end of simulation
							N_Start = 21
							N_End = 30
							@fastmath @inbounds @simd for iT=1:Nit
							# @fastmath @inbounds @simd for iZ=1:N_iRoot
								@fastmath @inbounds @simd for iZ=N_Start:N_End
									θroot_Mean[iOpt_Count] += (θ[iT,iZ] / hydro.θs[iZ]) * discret.ΔZ[iZ]
								end
							end
							# θroot_Mean[iOpt_Count] = θroot_Mean[iOpt_Count] / (Nit * sum(discret.ΔZ[1:N_iRoot]))
							θroot_Mean[iOpt_Count] = θroot_Mean[iOpt_Count] / (Nit * sum(discret.ΔZ[N_Start:N_End]))

						# Timing 
							Time_End = now()
							ΔRunTimeHypix[iOpt_Count] = value(Time_End - Time_Start) / 1000

						# Convergence rate
							iNonConverge_iOpt[iOpt_Count]          = iNonConverge
							Efficiency[iOpt_Count]                 = ceil(Int, cst.Day_2_Second * IterCount / ∑T[Nit] )
							ΔT_Average[iOpt_Count]                 = ceil(Int, mean(ΔT[i∑T_CalibrStart:Nit]))
							Global_WaterBalance[iOpt_Count]        = ∑∑WaterBalance
							Global_WaterBalance_NormPr[iOpt_Count] = 100.0 * ∑WaterBalance_η[Nit]
							∑∑ΔSink[iOpt_Count]                    = ∑ΔSink[Nit]
						
						# Ground water recharge
							∑ΔQ_Bot[iOpt_Count] = 0.0::Float64
							for iT=i∑T_CalibrStart:Nit
								∑ΔQ_Bot[iOpt_Count] = ∑ΔQ_Bot[iOpt_Count] + ΔT[iT] * Q[iT, Nz+1]
							end

							∑Q_Z = θaver.∑QZ(∑Q_Z, iOpt_Count, Nz, Z, ΔQ_Reduced; Zq=600.0)
				
						println("		=== ===START: summary  $iMultistep steps ...")

						if optionHypix.TopBoundary⍰ ≠ "Ψ"
							println("			∑PrThroughfall 			= ", round(∑Pr_Clim[iOpt_Count], digits=0), "  [mm]")
							println("			∑Pr_Soil 		= ", round(∑PrThroughfall[Nit] - ∑PrThroughfall[i∑T_CalibrStart], digits=0),  "  [mm]")
							println("			∑Pr_Intercepted/∑PrThroughfall 	= ", round(100. * (∑Pr_Clim[iOpt_Count] - (∑PrThroughfall[Nit]-∑PrThroughfall[i∑T_CalibrStart])) / (∑Pr_Clim[iOpt_Count] + eps()), digits=0),  "  [%]")
						end
						if optionHypix.RootWaterUptake
							println("			∑Pet_Net 		= ", ceil(Int, ∑Pet_Net[iOpt_Count]), "  [mm]")
							println("			∑Pet 			= ", ceil(Int, ∑Pet[Nit]- ∑Pet[i∑T_CalibrStart]), "  [mm]")
							println("			∑ΔSink/∑Pet_Net 	= ", ceil(Int, 100.0 * ∑ΔSink[Nit] /(∑Pet_Net[iOpt_Count] + eps(10.0))), "  [%]")
						end
						
						println("			∑SoilWaterContentRootEnd = ", round(θroot_Mean[iOpt_Count], digits=3), "  [mm]")
						println("			∑ΔSink 			= ", -ceil(Int, ∑∑ΔSink[iOpt_Count]), "  [mm]")
						println("			∑Infilt_Bot 		= ", -round(∑ΔQ_Bot[iOpt_Count],  digits=5), "  [mm]")
						println("			Hpond at end 		= ", ceil(Int, Hpond[Nit]), "  [mm] ")
						println("			∑Runof 				= ", round(∑Runoff_Reduced[end], digits=0), "  [mm]" )
		
						println("\n			Number_of_cells 	        = ", Nz, "[-]")
						println("			Global_WaterBalance_NormPr 	= ", round(Global_WaterBalance_NormPr[iOpt_Count], digits=8), "  [%]")
						println("			Global_WaterBalance 		= ", 	round(∑∑WaterBalance, digits=8), "  [mm]")
						println("			Average ΔT 			= ",  ΔT_Average[iOpt_Count] , "  [seconds]")
						println("			ΔTmin 				= ",   round(minimum(ΔT[i∑T_CalibrStart:Nit-1]), digits=0) , "  [seconds]")
						println("			ΔTmax 				= ",  round(maximum(ΔT[i∑T_CalibrStart:Nit-1]), digits=0) , "  [seconds]")
						println("			ΔT_HyPix 			= ", ceil(Int, ΔRunTimeHypix[iOpt_Count]) , "  [seconds]")

						println("\n			Efficiency 			= ", Efficiency[iOpt_Count], "  [iTer day-1]")
						println("			iNonConverge 			= ", iNonConverge, "  [count] \n")	

					# Computing average simulated θ to comapre it with average observed θ
					if optionHypix.θavr_RootZone && optionHypix.θobs	
						θsim_Aver = θaver.θAVER(discret; Z=Z, θ_Reduced=θ_Reduced, Nz=Nz, Nit_Reduced=Nit_Reduced, ZaverArray=[400.0])

						NoNaN =  .!( isnan.( θobs_Reduced[1:Nit_Reduced, 1]))

						θsim_Aver_Mean    = mean(θsim_Aver[NoNaN])
						θobs_Reduced_Mean = mean(filter(.!isnan, θobs_Reduced[NoNaN, 1]))

						θobs_Reduced_MeanAdj = θobs_Reduced[NoNaN, 1] .- θobs_Reduced_Mean
						θsim_Reduced_MeanAdj = θsim_Aver[NoNaN] .- θsim_Aver_Mean

						CccBest[iOpt_Count]      = stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(θobs_Reduced_MeanAdj, θsim_Reduced_MeanAdj)
						NseBest[iOpt_Count]      = stats.NSE(θobs_Reduced_MeanAdj, θsim_Reduced_MeanAdj)
						WilmotBest[iOpt_Count]   = stats.NSE_WILMOT(θobs_Reduced_MeanAdj, θsim_Reduced_MeanAdj)

					elseif  optionHypix.θobs
						for iZobs = 1:obsθ.Ndepth
                     CccBest[iOpt_Count]      += stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(θobs_Reduced[1:Nit_Reduced, iZobs], θ_Reduced[1:Nit_Reduced, obsθ.ithetaObs[iZobs]])
                     NseBest[iOpt_Count]      += stats.NSE(θobs_Reduced[1:Nit_Reduced, iZobs], θ_Reduced[1:Nit_Reduced, obsθ.ithetaObs[iZobs]])
                     WilmotBest[iOpt_Count]   += stats.NSE_WILMOT(θobs_Reduced[1:Nit_Reduced, iZobs], θ_Reduced[1:Nit_Reduced, obsθ.ithetaObs[iZobs]])
						end # loop

						θsim_Aver = Array{Float64, 2}(undef, 2, 2)
					else
						θsim_Aver = Array{Float64, 2}(undef, 2, 2)
					end

					if  optionHypix.θobs
						CccBest[iOpt_Count] = CccBest[iOpt_Count] / obsθ.Ndepth
						NseBest[iOpt_Count] = NseBest[iOpt_Count] / obsθ.Ndepth
						WilmotBest[iOpt_Count] = WilmotBest[iOpt_Count] / obsθ.Ndepth

						println("			Ccc 			= ", round(CccBest[iOpt_Count], digits=3))
						println("			Nse 			= ", round(NseBest[iOpt_Count], digits=3))
						println("			Wilmot                  = ", round(WilmotBest[iOpt_Count], digits=3))

						println("\n			Ccc_Average 	= ", round(mean(max.(CccBest[1:iScenario],0.0)), digits=3))
						println("			Nse_Average 	= ", round(mean(max.(NseBest[1:iScenario],0.0)), digits=3))
						println("			Wilmot_Average = ", round(mean(max.(WilmotBest[1:iScenario],0.0)), digits=3))
					end	
					println("		=== === END: summary \n")

					if optionHypix.Table
						if  optionHypix.θobs
							tableHypix.TABLE_HYPIX(∑∑ΔSink, ∑Pet_Net, ∑Pr_Clim, ∑Q_Z, ∑T, ∑WaterBalanceη_Reduced, ∑ΔQ_Bot, CccBest, Date_Reduced, discret, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, hydro, hydroHorizon, i∑T_CalibrStart, iMultistep, iNonConverge_iOpt, iOpt_Count, iScenario, N_SoilLayer, Nit, Nit_Reduced, NseBest, Nz, optionHypix, paramHypix, pathInputHypix, pathOutputHypix, SiteName, veg, WilmotBest, WofBest, Z, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPrGross_Reduced, ΔPrThroughfall_Reduced, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔRunTimeHypix, ΔSink_Reduced, ΔT_Average, θ_Reduced, θroot_Mean, Ψ_Reduced; θobs_Reduced=θobs_Reduced, θsim_Aver=θsim_Aver)
						else
							tableHypix.TABLE_HYPIX(∑∑ΔSink, ∑Pet_Net, ∑Pr_Clim, ∑Q_Z, ∑T, ∑WaterBalanceη_Reduced, ∑ΔQ_Bot, CccBest, Date_Reduced, discret, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, hydro, hydroHorizon, i∑T_CalibrStart, iMultistep, iNonConverge_iOpt, iOpt_Count, iScenario, N_SoilLayer, Nit, Nit_Reduced, NseBest, Nz, optionHypix, paramHypix, pathInputHypix, pathOutputHypix, SiteName, veg, WilmotBest, WofBest, Z, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPrGross_Reduced, ΔPrThroughfall_Reduced, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔRunTimeHypix, ΔSink_Reduced, ΔT_Average, θ_Reduced, θroot_Mean, Ψ_Reduced)
						end
					end
		
					if optionHypix.Ploting
					println("		=== === START: Plotting === ===")
				
					if optionHypix.θobs
						plotHypix.PLOT_HYPIX(∑T_Reduced, ∑ΔQ_Obs_Reduced, clim, Date_Reduced, dateHypix, discret, hydro, hydroHorizon, i∑T_CalibrStart_Day, iMultistep, iScenario, N_iRoot, N_Layer, Nit_Reduced, Nz, optionHypix, paramHypix, pathInputHypix, pathOutputHypix, SiteName, veg, Z, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPrGross_Reduced, ΔPrThroughfall_Reduced, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRootDensity, ΔRunoff_Reduced, ΔSink_Reduced, θ_Reduced; obsθ=obsθ, θobs_Reduced=θobs_Reduced, θsim_Aver=θsim_Aver)
					else
						plotHypix.PLOT_HYPIX(∑T_Reduced, ∑ΔQ_Obs_Reduced, clim, Date_Reduced, dateHypix, discret, hydro, hydroHorizon, i∑T_CalibrStart_Day, iMultistep, iScenario, N_iRoot, N_Layer, Nit_Reduced, Nz, optionHypix, paramHypix, pathInputHypix, pathOutputHypix, SiteName, veg, Z, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPrGross_Reduced, ΔPrThroughfall_Reduced, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRootDensity, ΔRunoff_Reduced, ΔSink_Reduced, θ_Reduced)
					end
		
					println("		=== === END: Plotting === === \n")
					end # if optionHypix.Plotting
			
			println("	=== === === END   ",iMultistep, "  steps ")
			println("	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n")
			end # for loop: iMultistep

		end # for iScenario = 1:N_Scenario
			
	end  # function: HYPIX_START
	# ------------------------------------------------------------------

end  # module hydro
# ............................................................

# @time hypixStart.HYPIX_START("LinkingFile.csv", "LYSIMETERS")
@time hypixStart.HYPIX_START("LinkingFile.csv", "ASHLEYOPT")
# hypixStart.HYPIX_START("LinkingFile.csv", "TESTCASE")
# @time hypixStart.HYPIX_START("LinkingFile_OVERSEER.csv", "SMAP")
# julia --check-bounds=no
