# =============================================================
#		MODULE: bestFunc
# =============================================================
module bestFunc
	import ..sorptivity, ..wrc, ..kunsat
	export  BEST_UNIVERSAL_START, CONVERT_3D_2_1D, CONVERT_1D_2_3D

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BEST_UNIVERSAL_START(∑Infilt_1D₀, ∑Infilt_3D₀, hydroInfilt₀, infiltParam₀, iZ, N_Infilt, option, Time₀; θini=infiltParam₀.θini[iZ])

			# Initializing
			Se_Ini = wrc.θ_2_Se(θ₁=θini, θs=hydroInfilt₀.θs[iZ], θr=hydroInfilt₀.θr[iZ])

			# Kr_θini = (kunsat.Se_2_KUNSAT(option.infilt, Se_Ini, iZ, hydroInfilt₀)) / hydroInfilt₀.Ks[iZ]
			Kr_θini = kunsat.KUNSAT_θΨSe(option.infilt, -1.0, iZ, hydroInfilt₀; Se₁=Se_Ini) / hydroInfilt₀.Ks[iZ]

			Sorptivity = sorptivity.SORPTIVITY(θini, iZ, hydroInfilt₀, option, option.infilt)

			A = bestFunc.A(θini, hydroInfilt₀.θs[iZ], iZ, infiltParam₀)

			B = bestFunc.B(iZ, Kr_θini, infiltParam₀)

			T_TransSteady = bestFunc.TIME_TRANS_STEADY(B, hydroInfilt₀.Ks[iZ], Sorptivity)

			for iT = 1:N_Infilt[iZ]
				if option.infilt.DataSingleDoubleRing⍰ == "Single"
					∑Infilt_3D₀[iZ, iT] = BEST_UNIVERSAL(iZ, A, B, Sorptivity, Time₀[iZ,iT], T_TransSteady, hydroInfilt₀, infiltParam₀, option)

				elseif option.infilt.DataSingleDoubleRing⍰ == "Double" 
					∑Infilt_1D₀[iZ, iT] = BEST_UNIVERSAL(iZ, A, B, Sorptivity, Time₀[iZ,iT], T_TransSteady, hydroInfilt₀, infiltParam₀, option) 
				end
			end  # for iT=1:N_Infilt[iZ] 

		return ∑Infilt_1D₀, ∑Infilt_3D₀, T_TransSteady
		end # function: BEST_UNIVERSAL_START

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST_UNIVERSAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BEST_UNIVERSAL(iZ, A, B, Sorptivity, Time₀, T_TransSteady, hydroInfilt₀, infiltParam₀, option)

			if option.infilt.DataSingleDoubleRing⍰ == "Single" #<>=<>=<>=<>=<>
				if Time₀ ≤ T_TransSteady
					return ∑Infilt_3D₀ = bestFunc.INFILTRATION_3D_TRANSIT(A, B, hydroInfilt₀.Ks[iZ], Sorptivity, Time₀)
				else
					return ∑Infilt_3D₀ = bestFunc.INFILTRATION_3D_STEADY(A, B, iZ, hydroInfilt₀.Ks[iZ], Sorptivity, Time₀, infiltParam₀, option, T_TransSteady)
				end # Time₀ <= T_TransSteady

			elseif option.infilt.DataSingleDoubleRing⍰ == "Double"  #<>=<>=<>=<>=<>
				if Time₀ ≤ T_TransSteady
					return ∑Infilt_1D₀ = bestFunc.INFILTRATION_1D_TRANSIT(B, hydroInfilt₀.Ks[iZ], Sorptivity, Time₀)
				else
					return ∑Infilt_1D₀ = bestFunc.INFILTRATION_1D_STEADY(B, iZ, hydroInfilt₀.Ks[iZ], Sorptivity, Time₀, infiltParam₀, option, T_TransSteady)
				end # Time₀ <= T_TransSteady
			end # option.∑Infilt_3D₀.Dimension	
		end  # function: BEST_UNIVERSAL
	#.................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_3D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_3D_TRANSIT(A, B, Ks, Sorptivity, Time₀)
			return Sorptivity * √ Time₀ + (A * (Sorptivity ^ 2.0) + B * Ks) * Time₀
		end  # function: INFILTRATION_3D_TRANSIT
	#.................................................................
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_1D_TRANSIT(B, Ks, Sorptivity, Time₀)
			return Sorptivity * √ Time₀ + B * Ks * Time₀
		end # function: INFILTRATION_1D_TRANSIT
	#.................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_3D_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_3D_STEADY(A, B, iZ, Ks, Sorptivity, Time₀, infiltParam₀, option, T_TransSteady)
			return bestFunc.INFILTRATION_3D_TRANSIT(A, B, Ks, Sorptivity, T_TransSteady) + (A * (Sorptivity ^ 2.0) + Ks) * (Time₀ - T_TransSteady)
		end  # function: INFILTRATION_3D_STEADY
	#.................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_1D_STEADY(B, iZ, Ks, Sorptivity, Time₀, infiltParam₀, option, T_TransSteady)
			return bestFunc.INFILTRATION_1D_TRANSIT(B, Ks, Sorptivity, T_TransSteady) + Ks * (Time₀ - T_TransSteady)
		end  # function: INFILTRATION_1D_STEADY
	#.................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERT_3D_2_1D
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CONVERT_3D_2_1D(∑Infilt_3D₀, ∑Infilt_1D₀, hydroInfilt₀, infiltParam₀, iZ, N_Infilt, option, Time₀; θini= infiltParam₀.θini[iZ])
				
			Sorptivity = sorptivity.SORPTIVITY(θini, iZ, hydroInfilt₀, option, option.infilt)

			A = bestFunc.A(θini, hydroInfilt₀.θs[iZ], iZ, infiltParam₀)

			for iT=1:N_Infilt[iZ]
				∑Infilt_1D₀[iZ,iT] = ∑Infilt_3D₀[iZ,iT] - A * Time₀[iZ,iT] * Sorptivity ^ 2.0
			end # iT

		return ∑Infilt_1D₀
		end  # function: CONVERT_3D_2_1D
	#.................................................................
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATIONtime_2_∑INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATIONtime_2_∑INFILT(hydroInfilt₀, infiltParam₀, iZ, option, Se_Ini; Infilt_10mm=10.0)

			# Initializing
				θini = wrc.Se_2_θ(Se₁=Se_Ini, θs=hydroInfilt₀.θs[iZ], θr=hydroInfilt₀.θr[iZ])

				Kr_θini = kunsat.KUNSAT_θΨSe(option.infilt, -1.0, iZ, hydroInfilt₀; Se₁=Se_Ini) / hydroInfilt₀.Ks[iZ]

				Sorptivity = sorptivity.SORPTIVITY(θini, iZ, hydroInfilt₀, option, option.infilt)

				A = bestFunc.A(θini, hydroInfilt₀.θs[iZ], iZ, infiltParam₀)

				B = bestFunc.B(iZ, Kr_θini, infiltParam₀)

				T_TransSteady = bestFunc.TIME_TRANS_STEADY(B, hydroInfilt₀.Ks[iZ], Sorptivity)

			# Iteration
         Time2Infiltrate = 0.0
         ∑Infilt_1D₂     = 0.0

				while ∑Infilt_1D₂ < Infilt_10mm
					Time2Infiltrate += 1.0 # Seconds

					if Time2Infiltrate ≤ T_TransSteady
						∑Infilt_1D₂ = bestFunc.INFILTRATION_1D_TRANSIT(B, hydroInfilt₀.Ks[iZ], Sorptivity, Time2Infiltrate)
					else
						∑Infilt_1D₂ = bestFunc.INFILTRATION_1D_STEADY(B, iZ, hydroInfilt₀.Ks[iZ], Sorptivity, Time2Infiltrate, infiltParam₀, option, T_TransSteady)
					end # Time₀ <= T_TransSteady
				end # iT

				# println("iZ=$iZ Se= $Se_Ini Time = $(Time2Infiltrate/60.0) minutes to reach Infilt_10mm= $Infilt_10mm")
	return Time2Infiltrate
	end  # function: name
	# ------------------------------------------------------------------
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERT_1D_2_3D
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CONVERT_1D_2_3D(∑Infilt_3D₀, ∑Infilt_1D₀, hydroInfilt₀, infiltParam₀, iZ, N_Infilt, option, Time₀; θini= infiltParam₀.θini[iZ])

			Sorptivity = sorptivity.SORPTIVITY(θini, iZ, hydroInfilt₀, option, option.infilt)

			A = bestFunc.A(θini, hydroInfilt₀.θs[iZ], iZ, infiltParam₀)

			for iT=1:N_Infilt[iZ]
				∑Infilt_3D₀[iZ,iT] = ∑Infilt_1D₀[iZ,iT] + A  * Time₀[iZ,iT] * Sorptivity ^ 2.0
			end # iT

			return ∑Infilt_3D₀
		end  # function: CONVERT_1D_2_3D
	#.................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_TRANS_STEADY 
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_TRANS_STEADY(B, Ks, Sorptivity)
			return (Sorptivity / (Ks * 2.0 * (1.0 - B)) ) ^ 2.0
		end # function: TIME_TRANS_STEADY
	#.................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : A
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function A(θini, θs, iZ, infiltParam₀)
			return  infiltParam₀.γ[iZ] / ( infiltParam₀.RingRadius[iZ] * (θs - θini)) # Units [mm-1]
		end  # function: A
	#.................................................................	


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : B
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function B(iZ, Kr_θini, infiltParam₀)
			return ((2.0 -  infiltParam₀.β[iZ]) + Kr_θini * (1.0 + infiltParam₀.β[iZ])) / 3.0
		end # function: B
	#.................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : C
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function C(B, infiltParam₀, iZ)
			return log(1.0 / infiltParam₀.β[iZ]) * (1.0 +  infiltParam₀.β[iZ]) / (6.0 * (1.0 -  infiltParam₀.β[iZ]) * (1.0 - B) )
		end # function: C
	#.................................................................

end # MODULE: bestFunc