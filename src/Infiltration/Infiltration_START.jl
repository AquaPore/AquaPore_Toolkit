# =============================================================
#		MODULE: infiltration
# =============================================================

include("OfBest.jl")
include("QuasiExact.jl")
include("TimeTransSteady.jl")
include("InfiltStruct.jl")
include("InfiltInitialize.jl")

module infiltStart
	import ..sorptivity, ..wrc, ..kunsat, ..infiltInitialize, ..bestFunc, ..stats, ..tool, ..quasiExact, ..ofBest, ..hydroRelation
	import BlackBoxOptim, Statistics
	export START_INFILTRATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function START_INFILTRATION(;∑Infilt_Obs, ∑Psd, hydro=[], hydroInfilt, infiltParam, N_Infilt, NiZ, option, param, Tinfilt)

		# INITIALIZE
		∑Infilt_1D,  ∑Infilt_3D, hydroInfilt, infiltOutput, Time = infiltInitialize.INFILT_INITIALIZE(∑Infilt_Obs, ∑Psd, hydroInfilt, infiltParam, N_Infilt, NiZ, option, param, Tinfilt)

		for iZ=1:NiZ
			println( "iZ= $iZ")

			# No optimization required running from hydro derived from laboratory
			if option.infilt.OptimizeRun⍰ == "Run" && option.run.HydroLabθΨ⍰ ≠ "No" #<>=<>=<>=<>=<>
				# Hydraulic param from laboratory
					hydroInfilt = deepcopy(hydro)
				 
				# Not to have errors
					infiltParam.θini[iZ] = max(hydroInfilt.θr[iZ] + eps(), infiltParam.θini[iZ])

				if option.infilt.Model⍰ == "Best_Univ"
					∑Infilt_1D, ∑Infilt_3D, T_TransStead =  bestFunc.BEST_UNIVERSAL_START(∑Infilt_1D, ∑Infilt_3D, hydroInfilt, infiltParam, iZ, N_Infilt, option, Time)

				elseif option.infilt.Model⍰ == "QuasiExact"
					∑Infilt_3D = quasiExact.HYDRO_2_INFILTRATION3D(∑Infilt_3D, hydroInfilt, infiltParam, iZ, N_Infilt, option, Time)

				end # option.infilt.Model⍰


			elseif option.infilt.OptimizeRun⍰ == "RunOptKs" && option.run.HydroLabθΨ⍰ ≠ "No" #<>=<>=<>=<>=<>	
				# Hydraulic param from laboratory
					hydroInfilt = deepcopy(hydro)
				
				# Not to have errors
				infiltParam.θini[iZ] = max(hydroInfilt.θr[iZ] + eps(), infiltParam.θini[iZ])
				
				SearchRange =[(log10(hydroInfilt.Ks_Min), log10(hydroInfilt.Ks_Max))]

				Optimization = BlackBoxOptim.bboptimize(P -> OF_INFILT_2_HYDRO(∑Infilt_1D, ∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, option, param, Time; Ks=10.0^P[1])[1]; SearchRange=SearchRange, NumDimensions=1, TraceMode=:silent)

				hydroInfilt.Ks[iZ] = 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[1]


			elseif option.infilt.OptimizeRun⍰ == "Opt" && option.infilt.HydroModel⍰ == "Kosugi" # <>=<>=<>=<>=<>	
				if option.infilt.σ_2_Ψm⍰ == "Constrained"
					SearchRange =[ (hydroInfilt.σ_Min[iZ], hydroInfilt.σ_Max[iZ]), (0.0, 1.0), (log10(hydroInfilt.Ks_Min[iZ]), log10(hydroInfilt.Ks_Max[iZ]))]

					Optimization = BlackBoxOptim.bboptimize(P -> OF_INFILT_2_HYDRO(∑Infilt_1D, ∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, option, param, Time; σ=P[1], Ψm=P[2], Ks=10.0^P[3])[1]; SearchRange=SearchRange, NumDimensions=3, TraceMode=:silent)

					hydroInfilt.σ[iZ]  = BlackBoxOptim.best_candidate(Optimization)[1]
					hydroInfilt.Ψm[iZ] = BlackBoxOptim.best_candidate(Optimization)[2]
					hydroInfilt.Ks[iZ] = 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[3]
		
					# Computing Ψm
					hydroInfilt = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydroInfilt, iZ, option.infilt, param.hydro)
						
				elseif option.infilt.σ_2_Ψm⍰ == "No"
					SearchRange =[ (hydroInfilt.σ_Min[iZ], hydroInfilt.σ_Max[iZ]), (log10(hydroInfilt.Ψm_Min[iZ]), log10(hydroInfilt.Ψm_Max[iZ])), (log10(hydroInfilt.Ks_Min[iZ]), log10(hydroInfilt.Ks_Max[iZ]))]

					Optimization = BlackBoxOptim.bboptimize(P -> OF_INFILT_2_HYDRO(∑Infilt_1D, ∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, option, param, Time; σ=P[1], Ψm=10.0^P[2], Ks=10.0^P[3])[1]; SearchRange=SearchRange, NumDimensions=3, TraceMode=:silent)

					hydroInfilt.σ[iZ]  = BlackBoxOptim.best_candidate(Optimization)[1]
					hydroInfilt.Ψm[iZ] = 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
					hydroInfilt.Ks[iZ] = 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[3]
				
				elseif option.infilt.σ_2_Ψm⍰ == "UniqueRelationship"
					SearchRange =[ (hydroInfilt.θr_Min[iZ], hydroInfilt.θr_Max[iZ]), (hydroInfilt.σ_Min[iZ], hydroInfilt.σ_Max[iZ]), (log10(hydroInfilt.Ks_Min[iZ]), log10(hydroInfilt.Ks_Max[iZ]))]

					Optimization = BlackBoxOptim.bboptimize(P -> OF_INFILT_2_HYDRO(∑Infilt_1D, ∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, option, param, Time; θr=P[1], σ=P[2], Ks=10.0^P[3])[1]; SearchRange=SearchRange, NumDimensions=3, TraceMode=:silent)

               hydroInfilt.θr[iZ] = BlackBoxOptim.best_candidate(Optimization)[1]
               hydroInfilt.σ[iZ]  = BlackBoxOptim.best_candidate(Optimization)[2]
               hydroInfilt.Ks[iZ] = 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[3]

					# Computing Ψm
					hydroInfilt = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydroInfilt, iZ, option.infilt, param.hydro)

				else
					error("option.infilt.σ_2_Ψm⍰ = $(option.infilt.σ_2_Ψm⍰) does not exist")
				
				end # option.infilt.σ_2_Ψm⍰
	
			else
				error("ERROR SoilWaterToolBox = $(option.infilt.Model⍰) not found")
			end # option.infilt.OptimizeRun⍰


			# OUTPUTS RUNNING THE OPTIMAL INFILTRATION
				infiltOutput.Sorptivity[iZ] = sorptivity.SORPTIVITY(infiltParam.θini[iZ], iZ, hydroInfilt, option, option.infilt) 

				if option.infilt.Model⍰ == "Best_Univ"
					∑Infilt_1D, ∑Infilt_3D, T_TransStead = bestFunc.BEST_UNIVERSAL_START(∑Infilt_1D, ∑Infilt_3D, hydroInfilt, infiltParam, iZ, N_Infilt, option, Time)

				elseif option.infilt.Model⍰ == "QuasiExact"
					∑Infilt_3D = quasiExact.HYDRO_2_INFILTRATION3D(∑Infilt_3D, hydroInfilt, infiltParam, iZ, N_Infilt, option, Time)
				end # option.infilt.Model⍰

			# Statistics
				iT_TransSteady = infiltOutput.iT_TransSteady_Data[iZ]

				infiltOutput.Nse_Trans[iZ] = 1.0 - stats.NSE_MINIMIZE( ∑Infilt_Obs[iZ,2:iT_TransSteady], ∑Infilt_3D[iZ, 2:iT_TransSteady]; Power=2.0)

				infiltOutput.Nse_Steady[iZ] = 1.0 - stats.NSE_MINIMIZE( log10.(∑Infilt_Obs[iZ,iT_TransSteady+1:N_Infilt[iZ]]), log10.(∑Infilt_3D[iZ,iT_TransSteady+1:N_Infilt[iZ]]); Power=2.0)

				infiltOutput.Nse[iZ] = 0.5 * infiltOutput.Nse_Trans[iZ] + 0.5 * infiltOutput.Nse_Steady[iZ]
		end # for iZ=1:NiZ


		# AVERAGE STATISTICS
         Nse_Trans  = Statistics.mean(infiltOutput.Nse_Trans[1:NiZ])
         Nse_Steady = Statistics.mean(infiltOutput.Nse_Steady[1:NiZ])
         Nse        = Statistics.mean(infiltOutput.Nse[1:NiZ])

			println("    ~ $(option.infilt.SorptivityModel⍰) ~")
			println("    ~ Infiltration model Nse= $Nse Nse_Trans= $Nse_Trans,  Nse_Steady= $Nse_Steady ~")


		# DERIVING INFILTRATION CURVES FOR DIFFERENT INITIAL Se & TIME TO INFILTRATE
		∑Infilt_1D_SeIni, Infilt_SeIni, Time_2_Infilt = INFILTRATION_SEini(hydroInfilt, infiltParam, N_Infilt, NiZ, option, Time)

		# ∑Infilt_1D_SeIni=[]
		# Infilt_SeIni = []
		# Time_2_Infilt =[10.0, 20.0]
		
		# CONVERTING INFILTRATION DIMENSIONS	 
			for iZ=1:NiZ
				if option.infilt.Model⍰ == "Best_Univ" && option.infilt.DataSingleDoubleRing⍰ == "Single"
					∑Infilt_1D = bestFunc.CONVERT_3D_2_1D(∑Infilt_3D, ∑Infilt_1D, hydroInfilt, infiltParam, iZ, N_Infilt, option, Time)
				
				elseif option.infilt.Model⍰ == "Best_Univ" && option.infilt.DataSingleDoubleRing⍰ == "Double"
					∑Infilt_3D = bestFunc.CONVERT_1D_2_3D(∑Infilt_3D, ∑Infilt_1D, hydroInfilt, infiltParam, iZ, N_Infilt, option, Time; θini= infiltParam.θini[iZ])

				elseif option.infilt.Model⍰ == "QuasiExact"  && option.infilt.DataSingleDoubleRing⍰ == "Single" 
					∑Infilt_1D = quasiExact.CONVERT_3D_2_1D(∑Infilt_3D, ∑Infilt_1D, hydroInfilt, infiltParam, iZ, N_Infilt, option, Time)

				else
					error("Not yet implemented option.infilt.Model⍰ == $(option.infilt.Model⍰)   && option.infilt.DataSingleDoubleRing⍰ == $(option.infilt.DataSingleDoubleRing⍰)")
				end
			end # for iZ=1:NiZ

	return ∑Infilt_1D, ∑Infilt_1D_SeIni, ∑Infilt_3D, hydroInfilt, Infilt_SeIni, infiltOutput, Time_2_Infilt
	end # FUNCTION: START_INFILTRATION


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_INFILT_2_HYDRO
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_INFILT_2_HYDRO(∑Infilt_1D, ∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, option, param, Time; σ=hydroInfilt.σ[iZ], Ψm=hydroInfilt.Ψm[iZ], θr=hydroInfilt.θr[iZ], θs=hydroInfilt.θs[iZ], Ks=hydroInfilt.Ks[iZ])

         hydroInfilt.θs[iZ]       = θs
         hydroInfilt.θr[iZ]       = θr
         hydroInfilt.Ks[iZ]       = Ks
         hydroInfilt.σ[iZ]        = σ
         hydroInfilt.Ψm[iZ]       = Ψm
         hydroInfilt.θsMacMat[iZ] = θs
         hydroInfilt.ΨmMac[iZ]    = Ψm
         hydroInfilt.σMac[iZ]     = σ

			hydroInfilt = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydroInfilt, iZ, option.infilt, param.hydro)

			if option.infilt.Model⍰ == "Best_Univ"
				return Nse = ofBest.OF_BEST(∑Infilt_1D, ∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, option, Time; W=0.5)

			elseif option.infilt.Model⍰ == "QuasiExact"
				return Nse =quasiExact.OF_QUASIEXACT(∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, option, Time)
				
			end # Option.infilt.Model

		end # FUNCTION: OF_INFILT_2_HYDRO_Best
	#.................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_SEINI
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_SEini(hydroInfilt₁, infiltParam, N_Infilt, NiZ, option, T₁; Infilt_SeIni = [0.0, 0.2, 0.4, 0.6, 0.8, 0.99])

			infiltParam₁ = deepcopy(infiltParam)

			# Initializing 
            N_Infilt_Max     = maximum(N_Infilt[1:NiZ])
            ∑Infilt_1D_SeIni = zeros(Float64, (NiZ, N_Infilt_Max, length(Infilt_SeIni)))
            Time_2_Infilt    = zeros(Float64, (NiZ, length(Infilt_SeIni)))
            ∑Infilt_3D₁      = zeros(Float64, (NiZ, N_Infilt_Max))
            ∑Infilt_1D₁      = zeros(Float64, (NiZ, N_Infilt_Max))

			# INFILTRATION 1D FOR DIFFERENT INITIAL SE

			for iZ=1:NiZ
				for (iSeIni, iiSeIni) in enumerate(Infilt_SeIni) 

					# Computing initial θini
						infiltParam₁.θini[iZ] = wrc.Se_2_θ(Se₁=iiSeIni, θs=hydroInfilt₁.θs[iZ], θr=hydroInfilt₁.θr[iZ])

					# Computing sorptivity
						# infiltOutput₁.Sorptivity[iZ] = sorptivity.SORPTIVITY(infiltParam₁.θini[iZ], iZ, hydroInfilt₁, option, option.infilt) 

					if option.infilt.Model⍰ == "Best_Univ"
						∑Infilt_1D₁, ∑Infilt_3D₁, ~ = bestFunc.BEST_UNIVERSAL_START(∑Infilt_1D₁, ∑Infilt_3D₁, hydroInfilt₁, infiltParam₁, iZ, N_Infilt, option, T₁)

						if option.infilt.DataSingleDoubleRing⍰ == "Single"
							# Converting to 1D
							∑Infilt_1D₁ = bestFunc.CONVERT_3D_2_1D(∑Infilt_3D₁, ∑Infilt_1D₁, hydroInfilt₁, infiltParam₁, iZ, N_Infilt, option, T₁)	
						end
						
					elseif option.infilt.Model⍰ == "QuasiExact"
						∑Infilt_3D₁ = quasiExact.HYDRO_2_INFILTRATION3D(∑Infilt_3D₁, hydroInfilt₁, infiltParam₁, iZ, N_Infilt, option, T₁)

						# Converting to 1D
						∑Infilt_1D₁ = quasiExact.CONVERT_3D_2_1D(∑Infilt_3D₁, ∑Infilt_1D₁, hydroInfilt₁, infiltParam₁, iZ, N_Infilt, option, T₁)
					end # option.infilt.Model⍰

				# Time to achieve 10mm of 1D ∑infiltration
					Time_2_Infilt[iZ, iSeIni] = bestFunc.INFILTRATIONtime_2_∑INFILT(hydroInfilt₁, infiltParam₁, iZ, option, iiSeIni; Infilt_10mm=10.0)

				# Saving
					∑Infilt_1D_SeIni[iZ, :, iSeIni] = ∑Infilt_1D₁[iZ , :]

				end # for iiSeIni_Vector in SeIni_Vector
			end # for iZ=1:NiZ

			infiltParam₁ = nothing
			
		return ∑Infilt_1D_SeIni, Infilt_SeIni, Time_2_Infilt
		end  # function: INFILTRATION_SEINI
	# ------------------------------------------------------------------

end # MODULE: infilt
	