# =============================================================
#		module: hypixOpt
# =============================================================
module hydrolabOpt
	import ..ofHydrolab, ..tool, ..optimize, ..hydroRelation, ..psdThetar, ..stats
	using BlackBoxOptim, Statistics
	export HYDROLABOPT_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIXOPT_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYDROLABOPT_START(;NiZ, ∑Psd, θ_θΨobs, Ψ_θΨobs, N_θΨobs, K_KΨobs=[0], Ψ_KΨobs=[0], N_KΨobs=1, hydro, hydroOther, option, optionₘ, optim, param, θϵ=0.005)
		for iZ = 1:NiZ
			# CORRECTION OF THE FEASIBLE RANGE ~~~
				θobs_Min = minimum(θ_θΨobs[iZ, 1:N_θΨobs[iZ]])  	# Smallest measure θ

				θobs_Max = maximum(θ_θΨobs[iZ, 1:N_θΨobs[iZ]])  	# Greatest measure θ

			# CORRECTING Θr ~~~
			# We overwrite  if optionₘ.θrOpt⍰=="ParamPsd"
				if ("θr" ∈ optim.ParamOpt)
					hydro.θr_Max[iZ] = max( min(θobs_Min-θϵ, hydro.θr_Max[iZ]), hydro.θr_Min[iZ] ) # Maximum value of θr

					# Changing the feasible range of θr
					iθr = findfirst(isequal("θr"), optim.ParamOpt)[1]
					optim.ParamOpt_Max[iθr] = hydro.θr_Max[iZ]

				elseif ("θr" ∉ optim.ParamOpt) && (optionₘ.θrOpt⍰=="ParamPsd") && (option.data.Psd) # Derive θr frpm PSD
					hydro.θr[iZ] = min(psdThetar.PSD_2_θr_FUNC(∑Psd, hydro, iZ, param), max(θobs_Min-θϵ, 0.0))
				end # if ("θr" ∈ optim.ParamOpt)

			# TEST IF EXIST Ψ=0  ~~~
				if minimum(Ψ_θΨobs[iZ,1:N_θΨobs[iZ]]) < eps(1000.0)
					🎏_0sOpt = false
				else
					🎏_0sOpt = true
				end

			# CORRECTING θS  ~~~
				if ("θs" ∈ optim.ParamOpt) && !🎏_0sOpt
					hydro.θs_Min[iZ] = θobs_Max * 0.75
					hydro.θs_Max[iZ] = θobs_Max * 1.1
					hydro.Φ[iZ] = θobs_Max / param.hydro.Coeff_Φ_2_θs

					# Changing the feasible range of θs
						iθs = findfirst(isequal("θs"), optim.ParamOpt)[1]
						optim.ParamOpt_Min[iθs] = hydro.θs_Min[iZ]
						optim.ParamOpt_Max[iθs] = hydro.θs_Max[iZ]

				elseif ("θs" ∉ optim.ParamOpt) && !🎏_0sOpt # <>=<>=<>=<>=<>
						hydro.θs[iZ] = θobs_Max
						hydro.Φ[iZ] = hydro.θs[iZ] / param.hydro.Coeff_Φ_2_θs
						
						# Reinforcing θs > θr
							hydro.θs[iZ] = max(hydro.θs[iZ], hydro.θr_Max[iZ] + 0.01)

				elseif  ("θs" ∉ optim.ParamOpt) # <>=<>=<>=<>=<>
					if hydro.Φ[iZ] *  param.hydro.Coeff_Φ_2_θs > θobs_Max + θϵ
						hydro.θs[iZ] = hydro.Φ[iZ] *  param.hydro.Coeff_Φ_2_θs
					elseif hydro.Φ[iZ] *  (param.hydro.Coeff_Φ_2_θs + 0.015) > θobs_Max + θϵ
						hydro.θs[iZ] = hydro.Φ[iZ] *  (param.hydro.Coeff_Φ_2_θs + 0.015)
					else
						hydro.θs[iZ] = max(hydro.Φ[iZ] - θϵ, θobs_Max + θϵ)
					end # hydro.Φ[iZ] * 0.95 > θobs_Max + θϵ

					# Reinforcing θs > θr
						hydro.θs[iZ] = max(hydro.θs[iZ], hydro.θr_Max[iZ] + 0.01)
				end
				
			# CORRECTING Ks  ~~~
				if option.data.Kθ
					if minimum(Ψ_KΨobs[iZ,1:N_KΨobs[iZ]]) < eps(100.0)
						🎏_KsOpt = false
					else
						🎏_KsOpt = true
					end
				end # if option.data.Kθ

				if option.data.Kθ && "Ks" ∈ optim.ParamOpt
					# test if exist Ψ=0
					if "Ks" ∈ optim.ParamOpt
						if minimum(Ψ_KΨobs[iZ,1:N_KΨobs[iZ]]) < eps(100.0)
							🎏_KsOpt = false
						else
							🎏_KsOpt = true
						end
					end # if "Ks" ∈ optim.ParamOpt

					K_KΨobs_Max = maximum(K_KΨobs[iZ, 1:N_KΨobs[iZ]])

					if 🎏_KsOpt && ("Ks" ∈ optim.ParamOpt)
						hydro.Ks_Min[iZ] = K_KΨobs_Max # Greatest measure of Kunsat)

						# Modifying the searchrange
						iKs = findfirst(isequal("Ks"), optim.ParamOpt)[1]
						optim.ParamOpt_Min[iKs] = hydro.Ks_Min[iZ]
						optim.ParamOpt_Max[iKs] = max(optim.ParamOpt_Max[iKs], hydro.Ks_Min[iZ] + 0.01)

					elseif !(🎏_KsOpt) && ("Ks" ∈ optim.ParamOpt)
						hydro.Ks_Max[iZ] = K_KΨobs_Max # Greatest measure of Kunsat

						# Modifying the searchrange
                     iKs                     = findfirst(isequal("Ks"), optim.ParamOpt)[1]
                     optim.ParamOpt_Max[iKs] = hydro.Ks_Max[iZ]
                     optim.ParamOpt_Min[iKs] = max(hydro.Ks_Max[iZ] - eps(1000.0), eps(10.0))

					elseif ("Ks" ∉ optim.ParamOpt)
                  hydro.Ks_Max[iZ] = K_KΨobs_Max
						hydro.Ks_Min[iZ] = K_KΨobs_Max
                  hydro.Ks[iZ]     = hydro.Ks_Max[iZ]

					end # "Ks" ∈ optim.ParamOpt
				end # if "Ks" ∈ optim.ParamOpt
			
			# Updated searchrange
				SearchRange = optimize.SEARCHRANGE(optionₘ, optim)


			# OPTIMIZATION: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				Optimization = BlackBoxOptim.bboptimize(X -> hydrolabOpt.OF_HYDROLAB(hydro, iZ, K_KΨobs, N_KΨobs, N_θΨobs, optim, option, optionₘ, param, X, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs); SearchRange=SearchRange, NumDimensions=optim.NparamOpt, TraceMode=:silent)

				X = BlackBoxOptim.best_candidate(Optimization)

				hydro = hydrolabOpt.PARAM_2_hydro(hydro, iZ, optim, optionₘ, param, X)

				# STATISTICS
					if option.data.Kθ && "Ks" ∈ optim.ParamOpt
						Of, Of_θΨ, Of_Kunsat = ofHydrolab.OF_WRC_KUNSAT(hydro, iZ, N_θΨobs, optim, optionₘ, θ_θΨobs, Ψ_θΨobs; K_KΨobs=K_KΨobs, N_KΨobs=N_KΨobs, Ψ_KΨobs=Ψ_KΨobs)
					else
						Of, Of_θΨ, Of_Kunsat = ofHydrolab.OF_WRC_KUNSAT(hydro, iZ, N_θΨobs, optim, optionₘ, θ_θΨobs, Ψ_θΨobs)
					end  

					hydroOther.Rmse[iZ], hydroOther.Rmse_KΨ[iZ], hydroOther.Rmse_θΨ[iZ] = ofHydrolab.OF_RMSE(option, optionₘ, iZ, θ_θΨobs, Ψ_θΨobs, N_θΨobs, K_KΨobs, Ψ_KΨobs, N_KΨobs, hydro, optim) 
		end # for iZ = 1:NiZ

		hydroOther.Nse_θΨ, ~, ~ = stats.NSE_θΨ(hydro, N_θΨobs, NiZ,  optionₘ, θ_θΨobs, Ψ_θΨobs)

		hydroOther.NseWilmot_θΨ, ~, ~ = stats.NSE_WILMOT_θΨ(hydro, N_θΨobs, NiZ,  optionₘ, θ_θΨobs, Ψ_θΨobs)
	
		if "Ks" ∈ optim.ParamOpt
			hydroOther.Nse_KΨ, ~, ~ = stats.NSE_KΨ(hydro, N_KΨobs, NiZ, optionₘ, K_KΨobs, Ψ_KΨobs)

			hydroOther.NseWilmot_KΨ, ~, ~ = stats.NSE_WILMOT_KΨ(hydro, N_KΨobs, NiZ, optionₘ, K_KΨobs, Ψ_KΨobs)

			hydroOther.Nse = (hydroOther.Nse_KΨ .+ hydroOther.Nse_θΨ) ./ 2.0
		else
			hydroOther.Nse = deepcopy(hydroOther.Nse_θΨ)
		end

		# OVERALL STATISTICS OF THE OPTIMIZATION
			Nse_θΨ_Aver = Statistics.mean(hydroOther.Nse_θΨ[1:NiZ])
			Nse_KΨ_Aver = Statistics.mean(max.(hydroOther.Nse_KΨ[1:NiZ], 0.0))

			NseWilmot_θΨ_Aver = Statistics.mean(hydroOther.NseWilmot_θΨ[1:NiZ])
			NseWilmot_KΨ_Aver = Statistics.mean(max.(hydroOther.NseWilmot_KΨ[1:NiZ], 0.0))

			Rmse_Aver    = Statistics.mean(hydroOther.Rmse[1:NiZ])
			Rmse_θΨ_Aver = Statistics.mean(hydroOther.Rmse_θΨ[1:NiZ])
			Rmse_KΨ_Aver = Statistics.mean(hydroOther.Rmse_KΨ[1:NiZ])
				
			if "Ks" ∈ optim.ParamOpt
				Nse_Aver = (Nse_θΨ_Aver + Nse_KΨ_Aver) / 2.0
			else
				Nse_Aver = Nse_θΨ_Aver
			end

			println("	=== === Optimizing Hydraulic parameters === ")
			println("    		~  Nse_θΨ= $(round(Nse_θΨ_Aver,digits=3)),  NseWilmot_θΨ= $(round(NseWilmot_θΨ_Aver,digits=3)), Nse_KΨ_Aver= $(round(Nse_KΨ_Aver,digits=3)), NseWilmot_KΨ= $(round(NseWilmot_KΨ_Aver,digits=3)), Nse = $(round(Nse_Aver,digits=3))  ~")
			println("    		~  Rmse_θΨ = $(round(Rmse_θΨ_Aver,digits=4)),  RmseLog10_KΨ = $(round(Rmse_KΨ_Aver,digits=4)), Rmse = $(round(Rmse_Aver,digits=4))  ~ \n")
			println( "	=== === ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ === ===")
	return hydro, hydroOther
	end  # function: HYPIXOPT_START


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_HYDROLAB(hydro, iZ, K_KΨobs, N_KΨobs, N_θΨobs, optim, option, optionₘ, param, X, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs)
			# New optimized which are put into the matching veg or hydro parameters
				hydro = hydrolabOpt.PARAM_2_hydro(hydro, iZ, optim, optionₘ, param, X)
		
			# Weighted Objective Function
			if option.data.Kθ && "Ks" ∈ optim.ParamOpt
				Of, Of_θΨ, Of_Kunsat = ofHydrolab.OF_WRC_KUNSAT(hydro, iZ, N_θΨobs, optim, optionₘ, θ_θΨobs, Ψ_θΨobs; K_KΨobs=K_KΨobs, N_KΨobs=N_KΨobs, Ψ_KΨobs=Ψ_KΨobs)
			else
				Of, Of_θΨ, Of_Kunsat = ofHydrolab.OF_WRC_KUNSAT(hydro, iZ, N_θΨobs, optim, optionₘ, θ_θΨobs, Ψ_θΨobs)
			end 
		return Of
		end  # function: OF_HYPIX


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PARAM_2_hydro(hydro, iZ, optim, optionₘ, param, X)
			for iParam = 1:optim.NparamOpt
				# Determening if parameters are Log transformed
					if (optim.ParamOpt_LogTransform[iParam]) && !(optim.ParamOpt[iParam]=="Ψm" && optionₘ.σ_2_Ψm⍰ == "Constrained")
						Paramₐ = expm1(X[iParam])
					else
						Paramₐ = X[iParam]
					end  # if: optim.ParamOpt_LogTransform

				# Getting the current values of every layer of the hydro parameter of interest
					vectParam = getfield(hydro, Symbol(optim.ParamOpt[iParam]))

				# Updating the value of the parameters for the soil wanting to optimize by keeping the values constant
					vectParam[iZ] = Paramₐ

				# Putting the updated hydro into hydro
					setfield!(hydro, Symbol(optim.ParamOpt[iParam]), vectParam)
			end # for loop

			# ==================== SPECIAL CASE ====================

			# RELATIONSHIP BETWEEN σ AND Ψm
			if (optionₘ.σ_2_Ψm⍰ ≠ "No") && ("Ψm" ∈ optim.ParamOpt)
				hydro = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydro, iZ, optionₘ, param.hydro; Pσ=3.0)
			elseif (optionₘ.σ_2_Ψm⍰ =="UniqueRelationship") 
				hydro = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydro, iZ, optionₘ, param.hydro; Pσ=3.0)

			end # optionₘ.σ_2_Ψm⍰ ≠ No

			#  <>=<>=<>=<>=<>=<> Relationship between σ and θr
			if optionₘ.θrOpt⍰=="σ_2_θr" && ("θr" ∉ optim.ParamOpt) && ("σ" ∈ optim.ParamOpt)
				hydro.θr[iZ] = hydroRelation.σ_2_θr(hydro, iZ)
			end

			# Converting θsMacMat_ƞ -> θsMacMat
			if  optionₘ.HydroModel⍰ == "Kosugi"
				hydro.θsMacMat[iZ] = hydro.θsMacMat_ƞ[iZ] * hydro.θs[iZ]
			end

		return hydro
		end  # function: PARAM

end  # module hypixOpt
# ............................................................