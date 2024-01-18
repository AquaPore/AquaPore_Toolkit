
# =============================================================
#		module: optIndivSoil
# =============================================================
module optIndivSoil
	import ..hydroRelation, ..ofHydrolab, ..optimize, ..psdThetar
	using BlackBoxOptim
	export OPTIMIZE_INDIVIDUALSOILS


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPTIMIZE_INDIVIDUALSOILS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function OPTIMIZE_INDIVIDUALSOILS(;∑Psd::Vector{Any}, hydro::Main.hydroStruct.KOSUGI, hydroOther::Main.hydroStruct.HYDRO_OTHER, K_KΨobs::Matrix{Float64}, N_KΨobs=1, N_θΨobs::Vector{Int64}, NiZ::Int64, Of_Sample::Vector{Float64}, optim::Main.reading.OPTIM, option::Main.options.OPTION, optionₘ::Main.options.HYDRO, param::Main.params.PARAM, θ_θΨobs::Matrix{Float64}, θϵ=0.005::Float64, Ψ_KΨobs::Matrix{Float64}, Ψ_θΨobs::Matrix{Float64})


		for iZ = 1:NiZ # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			# TEST IF EXIST Ψ ≈ 0  ~~~
				if minimum(Ψ_θΨobs[iZ, 1:N_θΨobs[iZ]]) < 0.01 
					🎏ΘΨ_0 = true
				else
					🎏ΘΨ_0 = false
				end

			# FEASIBLE RANGE OF θ_θΨobs ~~~
				θobs_Min = minimum(θ_θΨobs[iZ, 1:N_θΨobs[iZ]])  # Smallest measurement of θ
				θobs_Max = maximum(θ_θΨobs[iZ, 1:N_θΨobs[iZ]])  # Greatest measurement of θ


			# CORRECTING θr  ~~~~~
				hydro.θr_Max[iZ] = min(max(θobs_Min - θϵ, 0.0), hydro.θr_Max[iZ]) # Maximum value of θr

				if ("θr" ∈ optim.ParamOpt)
					# Changing the feasible range of θr
					iθr = findfirst(isequal("θr"), optim.ParamOpt)[1]
					optim.ParamOpt_Max[iθr] = hydro.θr_Max[iZ]
				end

				if ("θr" ∉ optim.ParamOpt) && (optionₘ.θrOpt⍰=="ParamPsd") && (option.data.Psd) # Derive θr frpm PSD
					hydro.θr[iZ] = psdThetar.PSD_2_θr_FUNC(∑Psd, hydro, iZ, param)
					hydro.θr[iZ] = max(min(hydro.θr[iZ], hydro.θr_Max[iZ]), hydro.θr_Min[iZ])
				end # if ("θr" ∈ optim.ParamOpt)


			# CORRECTING θs  ~~~~~
				if ("θs" ∈ optim.ParamOpt)
					if 🎏ΘΨ_0
						hydro.Φ[iZ] = θobs_Max / param.hydro.Coeff_Φ_2_θs

						hydro.θs_Min[iZ] = θobs_Max * 0.9
						hydro.θs_Max[iZ] = θobs_Max * 1.1

					elseif !(🎏ΘΨ_0) 
						hydro.θs_Min[iZ] = hydro.Φ[iZ] * 0.9
						hydro.θs_Max[iZ] = hydro.Φ[iZ] * 1.1
					end  # if: 🎏ΘΨ_0

					# Changing the feasible range of θs
						iθs = findfirst(isequal("θs"), optim.ParamOpt)[1]
						optim.ParamOpt_Min[iθs] = hydro.θs_Min[iZ]
						optim.ParamOpt_Max[iθs] = hydro.θs_Max[iZ]	

				elseif ("θs" ∉ optim.ParamOpt) # *****
					if 🎏ΘΨ_0
						hydro.θs[iZ] = θobs_Max
						hydro.Φ[iZ]  = hydro.θs[iZ] / param.hydro.Coeff_Φ_2_θs

					elseif !(🎏ΘΨ_0)
						hydro.θs[iZ] = max(hydro.Φ[iZ] * param.hydro.Coeff_Φ_2_θs, θobs_Max + 0.0015)

					end
				end # if "θs"


			# ~~~~ TEST IF EXIST K(Ψ=0)  ~~~
				if option.data.Kθ
					if minimum(Ψ_KΨobs[iZ,1:N_KΨobs[iZ]]) <  0.01 
						🎏_Ks = true
					else
						🎏_Ks = false
					end
					
					K_KΨobs_Max = maximum(K_KΨobs[iZ, 1:N_KΨobs[iZ]])
				end # if option.data.Kθ

				if option.data.Kθ && "Ks" ∈ optim.ParamOpt				
					if 🎏_Ks 
						hydro.Ks_Min[iZ] = K_KΨobs_Max * 0.95
						hydro.Ks_Max[iZ] = K_KΨobs_Max  * 1.1
					
					elseif !(🎏_Ks)
						hydro.Ks_Min[iZ] = max(hydro.Ks_Min[iZ], K_KΨobs_Max)
						hydro.Ks_Max[iZ] = max(hydro.Ks_Max[iZ], K_KΨobs_Max + 0.0015)
					
					end # "Ks" ∈ optim.ParamOpt

					# Modifying the searchrange
					iKs                     = findfirst(isequal("Ks"), optim.ParamOpt)[1]
					optim.ParamOpt_Min[iKs] = hydro.Ks_Min[iZ]
					optim.ParamOpt_Max[iKs] = hydro.Ks_Max[iZ]

				elseif option.data.Kθ && "Ks" ∉ optim.ParamOpt && 🎏_Ks
					hydro.Ks[iZ] = K_KΨobs_Max

				end # if option.data.Kθ && "Ks" ∈ optim.ParamOpt
			

			# OPTIMIZATION: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
				# Updated searchrange
				SearchRange = optimize.SEARCHRANGE(optionₘ, optim)

				Optimization = BlackBoxOptim.bboptimize(X -> optIndivSoil.OF_HYDROLAB(hydro, iZ, K_KΨobs, N_KΨobs, N_θΨobs, Of_Sample, optim, option, optionₘ, param, X, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs); SearchRange=SearchRange, NumDimensions=optim.NparamOpt, TraceMode=:silent)

				X = BlackBoxOptim.best_candidate(Optimization)

				hydro = optIndivSoil.PARAM_2_hydro(hydro, iZ, optim, optionₘ, param, X)

				# FINAL CORRECTION
					if optionₘ.σ_2_Ψm⍰ ≠ "No"
						hydro.ΨmacMat[iZ] = hydroRelation.FUNC_θsMacMatη_2_ΨmacMat(θs=hydro.θs[iZ], θsMacMat=hydro.θsMacMat[iZ], θr=hydro.θr[iZ], ΨmacMat_Max=hydro.ΨmacMat[iZ])
					end

				# STATISTICS
					if option.data.Kθ && "Ks" ∈ optim.ParamOpt
						Of_Sample, Of_θΨ, Of_Kunsat = ofHydrolab.OF_WRC_KUNSAT(hydro, iZ, N_θΨobs, Of_Sample, optim, optionₘ, θ_θΨobs, Ψ_θΨobs; K_KΨobs=K_KΨobs, N_KΨobs=N_KΨobs, Ψ_KΨobs=Ψ_KΨobs)
					else
						Of_Sample, Of_θΨ, Of_Kunsat = ofHydrolab.OF_WRC_KUNSAT(hydro, iZ, N_θΨobs, Of_Sample, optim, optionₘ, θ_θΨobs, Ψ_θΨobs)
					end  

					hydroOther.Rmse[iZ], hydroOther.Rmse_KΨ[iZ], hydroOther.Rmse_θΨ[iZ] = ofHydrolab.OF_RMSE(option, optionₘ, iZ, θ_θΨobs, Ψ_θΨobs, N_θΨobs, K_KΨobs, Ψ_KΨobs, N_KΨobs, hydro, optim)
		end # for iZ=1:NiZ
		
	return hydro, hydroOther, Of_Sample
	end  # function: OPTIMIZE_INDIVIDUALSOILS
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_HYDROLAB
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_HYDROLAB(hydro, iZ, K_KΨobs, N_KΨobs, N_θΨobs, Of_Sample, optim, option, optionₘ, param, X, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs)

			# New optimized which are put into the matching veg or hydro parameters
				hydro = optIndivSoil.PARAM_2_hydro(hydro, iZ, optim, optionₘ, param, X)
		
			# Weighted Objective Function
			if option.data.Kθ && "Ks" ∈ optim.ParamOpt
				Of_Sample, Of_θΨ, Of_Kunsat = ofHydrolab.OF_WRC_KUNSAT(hydro, iZ, N_θΨobs, Of_Sample, optim, optionₘ, θ_θΨobs, Ψ_θΨobs; K_KΨobs=K_KΨobs, N_KΨobs=N_KΨobs, Ψ_KΨobs=Ψ_KΨobs)
			else
				Of_Sample, Of_θΨ, Of_Kunsat = ofHydrolab.OF_WRC_KUNSAT(hydro, iZ, N_θΨobs, Of_Sample, optim, optionₘ, θ_θΨobs, Ψ_θΨobs)
			end 
		return Of_Sample[iZ]
		end  # function: OF_HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PARAM_2_hydro(hydro, iZ, optim, optionₘ, param, X; ΔMinΘₛ_Θᵣ=0.05)
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

			# ==================== SPECIAL CASES ====================

			# RELATIONSHIP BETWEEN ΨmacMat ➡ σMac & ΨmMac
				if optionₘ.ΨmacMat_2_σMac_ΨmMac
					ΨmacMat₁ = hydroRelation.FUNC_θsMacMatη_2_ΨmacMat( θs=hydro.θs[iZ], θsMacMat=hydro.θsMacMat[iZ], θr=hydro.θr[iZ], ΨmacMat_Max=hydro.ΨmacMat[iZ], ΨmacMat_Min=0.0, θsMacMat_η_Tresh=0.95) 

					hydro.σMac[iZ]  = hydroRelation.FUNC_ΨmacMat_2_σMac(ΨmacMat=ΨmacMat₁)
				
					hydro.ΨmMac[iZ] = hydroRelation.FUNC_ΨmacMat_2_ΨmMac(ΨmacMat=ΨmacMat₁, σMac=hydro.σMac[iZ])
				end

			# RELATIONSHIP BETWEEN σ AND Ψm
				if (optionₘ.σ_2_Ψm⍰ ≠ "No") && ("Ψm" ∈ optim.ParamOpt)
					hydro = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydro, iZ, optionₘ, param.hydro)

				elseif (optionₘ.σ_2_Ψm⍰ =="UniqueRelationship") 
					hydro = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydro, iZ, optionₘ, param.hydro)

				end # optionₘ.σ_2_Ψm⍰ ≠ No

			#  <>=<>=<>=<>=<>=<> Relationship between σ and θr
				if optionₘ.θrOpt⍰=="σ_2_θr" && ("θr" ∉ optim.ParamOpt) && ("σ" ∈ optim.ParamOpt)
					hydro.θr[iZ] = hydroRelation.σ_2_θr(hydro, iZ)
				end

			# Reinforcing θs >> Θr
				if hydro.θs[iZ] < hydro.θr[iZ] + ΔMinΘₛ_Θᵣ
					hydro.θr[iZ] = 0.0
				end

			# Converting θsMacMat_ƞ -> θsMacMat
				if  optionₘ.HydroModel⍰ == "Kosugi"
					hydro.θsMacMat[iZ] = min(hydro.θsMacMat_ƞ[iZ] * (hydro.θs[iZ] - hydro.θr[iZ]) + hydro.θr[iZ], hydro.θs[iZ])
				end

				if  optionₘ.HydroModel⍰ == "Kosugi"
					if hydro.θsMacMat[iZ] > hydro.θs[iZ]
						error("θsMacMat: $iZ $(hydro.θsMacMat[iZ])> $(hydro.θs[iZ])")
					end
				end

				if hydro.θr[iZ] > hydro.θs[iZ]
					error("Θr_θs: iZ = $iZ  , $(hydro.θr[iZ]) > $(hydro.θs[iZ])")
				end

				if hydro.θr[iZ] < 0.0
					error("Θr: $iZ $(hydro.θr[iZ]) < 0.0)")
				end

		return hydro
		end  # function: PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end # module optIndivSoil