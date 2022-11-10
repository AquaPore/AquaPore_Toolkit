module interception

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RAINFALL_INTERCEPTION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RAINFALL_INTERCEPTION_START(∑Pet_Climate::Vector{Float64}, ∑PrThroughfall_Climate::Vector{Float64}, clim, optionHypix, veg)
				
			# INTERCEPTION MODEL
				Sint               = 0.0::Float64
				clim.Pr_Through[1] = 0.0::Float64
				clim.Pet[1]        = 0.0::Float64

			@fastmath for iT = 2:clim.N_Climate
				# Maximum water storage of the vegetation
					if optionHypix.Lai_2_SintMax
						Sint_Sat = LAI_2_SINTMAX(clim, iT, veg)
					else
						Sint_Sat = veg.Sint_Sat
					end
					
				# GapFraction: fraction of precipitation that reaches the ground surface through gaps in the canopy
					GapFrac = LAI_2_GAPFRAC(clim, iT, veg)
					
				# Rainfall interception model
					Sint, ΔEvap_Int, ΔPr_Through = RAINFALL_INTERCEPTION(GapFrac, Sint, Sint_Sat, clim.Pet[iT], clim.Pr[iT])

					clim.Pr_Through[iT] = ΔPr_Through

					∑PrThroughfall_Climate[iT] = ∑PrThroughfall_Climate[iT-1] + ΔPr_Through

					∑Pet_Climate[iT] = ∑Pet_Climate[iT-1] + clim.Pet[iT] - ΔEvap_Int
			end  # for  iT = 1:clim.N_Climate 

		return ∑Pet_Climate, ∑PrThroughfall_Climate, clim
		end  # function: RAINFALL_INTERCEPTION
	#-----------------------------------------------------------------
	

	"""
		RAINFALL_INTERCEPTION

	ΔPr = Precipitation which falls on top of canopy [mm / Δtime] 
	ΔPet_Intercept =  Potential Evapotranspiration
	Sint = current storage of water in the vegetation [mm]
	Sint_Sat = Maximum storage of water in the vegetation [mm]
	GapFrac = % of ΔPr which is intercepted by the vegetation at the very beginning [0-1]
	"""
		function RAINFALL_INTERCEPTION(GapFrac::Float64, Sint::Float64, Sint_Sat::Float64, ΔPet_Int::Float64, ΔPr::Float64; PevapInt = 0.666)
			# Requires initial Sint

			# ΔPr: falls on top of the vegetation; ΔPr_Int: amount of ΔPr which is intercepted by the vegetaion which depends on parameter GapFrac; ΔPr_Ground: precipitation which is not intercepted
				ΔPr_Ground = GapFrac * ΔPr
		
				ΔPr_Int = ΔPr - ΔPr_Ground
			
			# New Sint 
				Sint₂ = min(Sint + ΔPr_Int, Sint_Sat)
			
			# According to Rutter et al. (1971), evaporation from wet canopies is assumed to be proportional to the fraction of the canopy that is wet (Sint / Sint_Sat) computed following Deardorff (1978):
				ΔEvap_Int = min(ΔPet_Int * (Sint₂ / Sint_Sat) ^ PevapInt, Sint₂)
		
			# ΔPr that overflows because the vegetation cannot store more water
				ΔPr_Over = max(Sint + ΔPr_Int - ΔEvap_Int - Sint_Sat, 0.0::Float64)

			# Amount of water stored in the vegetaion. Sint =< Sint_Sat
				Sint = Sint + ΔPr_Int - ΔEvap_Int - ΔPr_Over

			# Total amount of throughfall
				ΔPr_Through = ΔPr_Ground + ΔPr_Over

		return Sint, ΔEvap_Int, ΔPr_Through
		end # RAINFALL_INTERCEPTION
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : LAI_2_SINTMAX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		"""
		LAI_2_SINTMAX

		Converts LAI to Sint_Sat
		The general equation based on Menzel’s equation
		"""
		function LAI_2_SINTMAX(clim, iT::Int64, veg)

			return veg.Sint_Lai * log(1.0 + clim.Lai[iT])
		end # function: LAI_2_storageVeg_Max
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : LAI_2_GAPFRAC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function LAI_2_GAPFRAC(clim, iT::Int64, veg)
			return 1.0 - exp(-veg.ExtinctCoefRadiation * clim.Lai[iT])
			#  println("			~ GapFrac = ", round(GapFrac, digits=2),"  ~")
		end  # function: LAI_2_GAPFRAC
	#-----------------------------------------------------------------

end # module interception