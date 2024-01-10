# =============================================================
#		module: ΨminΨmax
# =============================================================
module ΨminΨmax

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΨMINΨMAX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ΨMINΨMAX(θs, θsMacMat, σ, σMac, Ψm, ΨmMac; Pσ₂ =5.0)
			# Ψ_Min
				# Have we got macropore ?
				if θs - θsMacMat > 0.001
					Ψ_Min = exp( log(ΨmMac) - σMac * Pσ₂) 
				else
					Ψ_Min = exp( log(Ψm) - σ * Pσ₂)
				end  # if: θs -θsMacMat > 0.01

			#  Ψ_Max
				Ψ_Max = exp(log(Ψm) + σ * Pσ₂)
		return Ψ_Max, Ψ_Min
		end  # function: ΨMINΨMAX
		
	end  # module ΨminΨmax
# ............................................................