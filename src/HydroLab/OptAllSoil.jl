# =============================================================
#		module: optAllSoil
# =============================================================
module optAllSoil
import ..ofHydrolab, ..optimize, ..optIndivSoil
using BlackBoxOptim
export OPTIMIZE_ALLSOILS

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPTIMIZE_ALLSOILS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function OPTIMIZE_ALLSOILS(;∑Psd, hydro, hydroOther, K_KΨobs=[0], N_KΨobs=1, N_θΨobs, NiZ, optim, optimAllSoils, option, optionₘ, param, θ_θΨobs, Ψ_KΨobs=[0], Ψ_θΨob)


		# OPTIMISATION HYDRAULIC PARAMETERS ALL SOILS INDIVIDUALLY AND ALL SOILS
		hydro, hydroOther, Of_Sample = optIndivSoil.OPTIMIZE_INDIVIDUALSOILS(;∑Psd, hydro, hydroOther, K_KΨobs, N_KΨobs, N_θΨobs, NiZ, Of_Sample, optim, option, optionₘ, param, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs)
		
	return
	end  # function: OPTIMIZE_ALLSOILS
	# ------------------------------------------------------------------
	
end  # module optAllSoil
# ............................................................