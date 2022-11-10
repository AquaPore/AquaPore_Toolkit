# =============================================================
#		MODULE: et
# =============================================================
module pet
	export BEER_LAMBERT_LAW

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  BEER_LAMBERT_LAW
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BEER_LAMBERT_LAW(iT, iT_Pr, ΔPet, veg, clim)

			ΔPet_Evap = ΔPet[iT] * exp(-clim.Lai[iT_Pr-1] * veg.ExtinctCoefRadiation)
			
			ΔPet_Transp = ΔPet[iT] - ΔPet_Evap

		return ΔPet_Evap, ΔPet_Transp
		end  # function BEER_LAMBERT_LAW
	#-----------------------------------------------------------------
	
end  # module et
# ............................................................