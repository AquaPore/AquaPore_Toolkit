module quasiExact 
	import ..sorptivity, ..wrc, ..kunsat
	import BlackBoxOptim, Optim
 	export CONVERT_3D_2_1D, HYDRO_2_INFILTRATION3D, OF_QUASIEXACT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_2_TIMEη
	#		= COMPUTE NORMALISED TIME: Tη =
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_2_TIMEη(Sorptivity, Time₀, ΔK)
			return Time_η = Time₀ * 2.0 * (ΔK / Sorptivity) ^ 2.0
		end # function: TIME_2_TIMEη


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION3D_2_1D
	# 		= TRANSFORMS INFILTRATION_3D TO INFILTRATION_1D =
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function CONVERT_3D_2_1D(∑Infilt_3D₀, ∑Infilt_1D₀, hydroInfilt₀, infiltParam₀, iZ, N_Infilt, option, Time₀; θini= infiltParam₀.θini[iZ])

		Δθ = hydroInfilt₀.θs[iZ] - θini

		Sorptivity = sorptivity.SORPTIVITY(θini, iZ, hydroInfilt₀, option, option.infilt)

		for iT = 1:N_Infilt[iZ]
			∑Infilt_1D₀[iZ,iT] = ∑Infilt_3D₀[iZ,iT] - (Time₀[iZ,iT] * infiltParam₀.γ[iZ] * Sorptivity ^ 2.0) / (infiltParam₀.RingRadius[iZ] * Δθ)
		end
	return ∑Infilt_1D₀
	end # function : INFILTRATION3D_2_1D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  INFILTTRATIONη_2_3D
	# 		TRANSFORMS NORMALIZED INFILTRATION TO INFILTRATION-3D
	#		Function compute infiltration-1d from normalized infiltration
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTTRATIONη_2_3D(Infilt_η, infiltParam₀, iZ, K_θini, Sorptivity, Time₀, Time_η, ΔK, Δθ)
			
			INFILTTRATIONη_2_1D(Infilt_η, K_θini, Sorptivity, Time₀, ΔK) = K_θini * Time₀ + Infilt_η * (Sorptivity ^ 2.0) / (2.0 * ΔK)

			ΔI_η = Time_η * infiltParam₀.γ[iZ]

		return ∑Infilt_3D₀ = INFILTTRATIONη_2_1D(Infilt_η, K_θini, Sorptivity, Time₀, ΔK) + ΔI_η * (Sorptivity ^ 4.0) / (2.0 * infiltParam₀.RingRadius[iZ] * Δθ * (ΔK ^ 2.0))
		end # function :  INFILTTRATIONη_2_3D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION3D_2_η
 	# 		TRANSFORMS INFILTRATION-3D TO NORMALIZED INFILTRATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION3D_2_η(∑Infilt_3D₀, infiltParam₀, iZ, K_θini, Sorptivity, Time₀, ΔK, Δθ; ϵ=eps())
			return Infilt_η = max((2.0 * ΔK / Sorptivity ^ 2.0) * (∑Infilt_3D₀ - K_θini * Time₀ - infiltParam₀.γ[iZ] * TIME_2_TIMEη(Sorptivity, Time₀, ΔK) * (Sorptivity^4.0) / (infiltParam₀.RingRadius[iZ] * Δθ * 2.0* (ΔK^2.0)) ), ϵ)
		end # function INFILTRATION3D_2_η


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDRO_2_INFILTRATION3D
	# 		COMPUTE INFILTRATION_3D FROM OPTIMIZED HYDRAULIC PARAMETERS
	# 		Solving quasiexact solution
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
		function HYDRO_2_INFILTRATION3D(∑Infilt_3D₀, hydroInfilt₀, infiltParam₀, iZ, N_Infilt, option, Time₀; Infilt_η_Max_Start=0.5)

			Infilt_η = fill(0.0::Float64, N_Infilt[iZ])

			Sorptivity = sorptivity.SORPTIVITY(infiltParam₀.θini[iZ], iZ, hydroInfilt₀, option, option.infilt)

			Se_Ini = wrc.θ_2_Se(θ₁=infiltParam₀.θini[iZ], θs=hydroInfilt₀.θs[iZ], θr=hydroInfilt₀.θr[iZ])

			# K_θini = kunsat.Se_2_KUNSAT(option.infilt, Se_Ini, iZ, hydroInfilt₀)
			K_θini = kunsat.KUNSAT_θΨSe(option.infilt, -1.0, iZ, hydroInfilt₀; θ₁=-1.0, Se₁=Se_Ini)

			ΔK = hydroInfilt₀.Ks[iZ] - K_θini

			Δθ = hydroInfilt₀.θs[iZ] - infiltParam₀.θini[iZ]
		
			# At t=1
            ∑Infilt_3D₀[1] = 0.0
            Infilt_η[1]   = 0.0
            Infilt_η_Min  = 10^-8
            Infilt_η_Max  = Infilt_η_Max_Start #Since Time₀[1] = 0

			# ~~~~~~~~~~~~~~~~~~~~
			function OF_QUASIEXACTη(Infilt_η, infiltParam₀, iZ, Time_η)
				Left_Term = Time_η

				Right_Term = (1.0 / (1.0 - infiltParam₀.β[iZ])) * (Infilt_η - log((exp(infiltParam₀.β[iZ] * Infilt_η) + infiltParam₀.β[iZ] - 1.0) / infiltParam₀.β[iZ]))
				
				if Right_Term < 0.0
					return OF = 10000.0 * exp(Infilt_η)
				else
					return OF = abs(Left_Term - Right_Term)
				end
			end # function OF_QUASIEXACTη ~~~~~~~~~~~~~~~~~~~~


			for iT = 2:N_Infilt[iZ] # Looping for every time step
				Time_η = TIME_2_TIMEη(Sorptivity, Time₀[iZ,iT], ΔK)

				# Solving for Infilt_η
					Optimization = Optim.optimize(Infilt_η -> OF_QUASIEXACTη(Infilt_η, infiltParam₀, iZ, Time_η), Infilt_η_Min, Infilt_η_Max, Optim.GoldenSection())

					Infilt_η[iT] = Optim.minimizer(Optimization)[1]

				# Deriving the new bounds such that infiltration increases with time & the slope decreases with time
					Infilt_η_Min = Infilt_η[iT]

				# Maximum infiltration rate for Time₀+1: (Infiltration[T2] - Infiltration[T1]) / (T2 - T1) which is 1 seconds
					if iT <= N_Infilt[iZ] - 1
						Infilt_η_Max = Infilt_η[iT] + (Time₀[iZ,iT+1]- Time₀[iZ,iT]) * (Infilt_η[iT] - Infilt_η[iT-1]) / (Time₀[iZ,iT] - Time₀[iZ,iT-1])
					else
						Infilt_η_Max = Infilt_η[iT] + (Infilt_η[iT] - Infilt_η[iT-1])
						
					end

				# Transforming INFILTTRATIONη to INFILTRATION3D 
					∑Infilt_3D₀[iZ,iT] =  INFILTTRATIONη_2_3D(Infilt_η[iT], infiltParam₀, iZ, K_θini, Sorptivity, Time₀[iZ,iT], Time_η, ΔK, Δθ)
			end
		return ∑Infilt_3D₀
		end # function: HYDRO_2_INFILTRATION3D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_QUASIEXACT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	function OF_QUASIEXACT(∑Infilt_Obs, hydroInfilt₀, infiltOutput₀, infiltParam₀, iZ, N_Infilt, option, Time₀; W=0.8)

		Se_Ini = wrc.θ_2_Se(θ₁=infiltParam₀.θini[iZ], θs=hydroInfilt₀.θs[iZ], θr=hydroInfilt₀.θr[iZ])
		
		# K_θini = kunsat.Se_2_KUNSAT(option.infilt, Se_Ini, iZ, hydroInfilt₀)
		K_θini = kunsat.KUNSAT_θΨSe(option.infilt, -1.0, iZ, hydroInfilt₀; θ₁=-1.0, Se₁=Se_Ini)

		Kr_θini = K_θini / hydroInfilt₀.Ks[iZ]
		
		ΔK = hydroInfilt₀.Ks[iZ] - K_θini
		
		Δθ = hydroInfilt₀.θs[iZ] - infiltParam₀.θini[iZ]
		
		Sorptivity = sorptivity.SORPTIVITY(infiltParam₀.θini[iZ], iZ, hydroInfilt₀, option, option.infilt)

		Left_Term = zeros(Float64, N_Infilt[iZ])

		Right_Term = zeros(Float64, N_Infilt[iZ])

		iT_TransSteady = infiltOutput₀.iT_TransSteady_Data[iZ]
		
		Of_Penalty = 0.0 ; Of_Stead = 0.0 ; Of_Trans = 0.0
		
		for iT = 2:N_Infilt[iZ]
			Time_η = TIME_2_TIMEη(Sorptivity, Time₀[iZ,iT], ΔK)

			Infilt_η = INFILTRATION3D_2_η(∑Infilt_Obs[iZ,iT], infiltParam₀, iZ, K_θini, Sorptivity, Time₀[iZ,iT], ΔK, Δθ)

			Left_Term[iT] = Time_η

			Right_Term[iT] = (1.0 / (1.0 - infiltParam₀.β[iZ])) * (Infilt_η - log((exp(infiltParam₀.β[iZ] * Infilt_η) + infiltParam₀.β[iZ] - 1.0) / infiltParam₀.β[iZ]))
			
			if Right_Term[iT] > 0.0
				if iT <= infiltOutput₀.iT_TransSteady_Data[iZ]
					Of_Trans += ((Left_Term[iT]) - (Right_Term[iT])) ^ 2.0
				else
					Of_Stead += (log10(Left_Term[iT]) - log10(Right_Term[iT])) ^ 2.0
				end #  iT <= infiltOutput₀.iT_TransSteady_Data
			else
				Of_Penalty += 1000.0 * exp(Infilt_η)
				Right_Term[iT] = 0.0
			end #  Right_Term[iT] > 0.0
		end #  Right_Term[iT] > 0.0

	return Wof = (W * Of_Trans / Float64(iT_TransSteady-1)) + ((1.0 - W) * Of_Stead / Float64(N_Infilt[iZ] - iT_TransSteady + 1)) + Of_Penalty
	end # function: OF_INFILT_2_HYDRO
	


end # module quasiExact