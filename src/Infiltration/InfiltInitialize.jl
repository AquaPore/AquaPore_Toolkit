# =============================================================
#		MODULE: HYDRO INITIALIZE
# =============================================================
module infiltInitialize
	import ..hydroStruct, ..psdThetar, ..timeTransSteady, ..infiltStruct
	export INFILT_INITIALIZE

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILT_INITIALIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function INFILT_INITIALIZE(∑Infilt_Obs, ∑Psd, hydroInfilt, infiltParam, N_Infilt, NiZ, option, param, Tinfilt)

		# Computing ∑Infilt_1D at different SeIni
		SeIni_Vector = [0.0, 0.25, 0.5, 0.75]

		# CORRECTION FOR θr & θs
			for iZ=1:NiZ
				# θr computation
					if option.run.IntergranularMixingPsd
						hydroInfilt.θr[iZ] = psdThetar.PSD_2_θr_FUNC(∑Psd, hydroInfilt, iZ, param)

					else
						hydroInfilt.θr[iZ] = 0.0

					end # option.run.IntergranularMixingPsd

				# θr < θini
					infiltParam.θini[iZ] = max(hydroInfilt.θr[iZ] + eps(), infiltParam.θini[iZ])

				# θs computation
					hydroInfilt.θs[iZ] = param.hydro.Coeff_Φ_2_θs * hydroInfilt.Φ[iZ]

				# Checking for θini
					infiltParam.θini[iZ] = min(infiltParam.θini[iZ], hydroInfilt.θs[iZ] - 0.0015) 
			end  # for iZ=1:NiZ

		# TIME FLUX CORRECTION
			N_Infilt_Max = maximum(N_Infilt[1:NiZ])

			T = fill(0.0::Float64, (NiZ, N_Infilt_Max))
			for iZ=1:NiZ
				# T[iZ,1] = ( (Tinfilt[iZ,1] ^0.5 + (0.0)^0.5) / 2.0 ) ^ 2.0
				T[iZ,1] = Tinfilt[iZ,1]

				for iInfilt=2:N_Infilt[iZ]
					T[iZ,iInfilt] = ( (Tinfilt[iZ,iInfilt] ^ 0.5 + Tinfilt[iZ,iInfilt-1] ^ 0.5) / 2.0 ) ^ 2.0
				end	
			end #  iZ=1:NiZ

		# STRUCTURE OF INFILTRATION PARAMETERS
			infiltOutput = infiltStruct.INFILTSTRUCT(NiZ)

		# DETERMENING WHEN STEADY STATE OCCURES
			infiltOutput = timeTransSteady.∑INFIlT_2_TIMETRANSSTEADY(∑Infilt_Obs, infiltOutput, N_Infilt, NiZ, param, T) 

		# Initializing Infiltration		
			∑Infilt_3D = fill(0.0::Float64, (NiZ, N_Infilt_Max))
			∑Infilt_1D = fill(0.0::Float64, (NiZ, N_Infilt_Max))
			∑Infilt_1D_SeIni = fill(0.0::Float64, (NiZ, N_Infilt_Max, length(SeIni_Vector)))

	return ∑Infilt_1D, ∑Infilt_1D_SeIni, ∑Infilt_3D, hydroInfilt, infiltOutput, SeIni_Vector, T
	end  # function: INFILT_INITIALIZE


end # module hydroInitialize
# ............................................................