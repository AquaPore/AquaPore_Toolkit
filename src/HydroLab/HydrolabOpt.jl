# =============================================================
#		module: hypixOpt
# =============================================================
module hydrolabOpt

	import   ..stats, ..optIndivSoil
	using  Statistics
	export HYDROLABOPT_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIXOPT_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYDROLABOPT_START(;∑Psd, hydro, hydroOther, K_KΨobs=[0], N_KΨobs=1, N_θΨobs, NiZ, optim, optimAllSoils, option, optionₘ, param, θ_θΨobs, Ψ_KΨobs=[0], Ψ_θΨobs)
		
		# Initiating arrays 
			Of_Sample = zeros(Float64, NiZ)

			if optimAllSoils.🎏_Opt
			# OPTIMISATION INDIVDUAL
				hydro, hydroOther, Of_Sample = optIndivSoil.OPTIMIZE_INDIVIDUALSOILS(;∑Psd, hydro, hydroOther, K_KΨobs, N_KΨobs, N_θΨobs, NiZ, Of_Sample, optim, option, optionₘ, param, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs)
		else
		# OPTIMISATION HYDRAULIC PARAMETERS ALL SOILS INDIVIDUALLY AND ALL SOILS
			hydro, hydroOther, Of_Sample = optIndivSoil.OPTIMIZE_INDIVIDUALSOILS(;∑Psd, hydro, hydroOther, K_KΨobs, N_KΨobs, N_θΨobs, NiZ, Of_Sample, optim, option, optionₘ, param, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs)
		end

		# STATISCS
         hydroOther.Nse_θΨ, ~, ~       = stats.NSE_θΨ(hydro, N_θΨobs, NiZ,  optionₘ, θ_θΨobs, Ψ_θΨobs)
         hydroOther.NseWilmot_θΨ, ~, ~ = stats.NSE_WILMOT_θΨ(hydro, N_θΨobs, NiZ,  optionₘ, θ_θΨobs, Ψ_θΨobs)
		
			if "Ks" ∈ optim.ParamOpt
            hydroOther.Nse_KΨ, ~, ~       = stats.NSE_KΨ(hydro, N_KΨobs, NiZ, optionₘ, K_KΨobs, Ψ_KΨobs)
            hydroOther.NseWilmot_KΨ, ~, ~ = stats.NSE_WILMOT_KΨ(hydro, N_KΨobs, NiZ, optionₘ, K_KΨobs, Ψ_KΨobs)
            hydroOther.Nse                = (hydroOther.Nse_KΨ .+ hydroOther.Nse_θΨ) ./ 2.0

			else
				hydroOther.Nse = deepcopy(hydroOther.Nse_θΨ)
			end

		# OVERALL STATISTICS OF THE OPTIMIZATION
         Nse_θΨ_Aver       = Statistics.mean(hydroOther.Nse_θΨ[1:NiZ])
         Nse_KΨ_Aver       = Statistics.mean(max.(hydroOther.Nse_KΨ[1:NiZ], 0.0))

         NseWilmot_θΨ_Aver = Statistics.mean(hydroOther.NseWilmot_θΨ[1:NiZ])
         NseWilmot_KΨ_Aver = Statistics.mean(max.(hydroOther.NseWilmot_KΨ[1:NiZ], 0.0))

         Rmse_Aver         = Statistics.mean(hydroOther.Rmse[1:NiZ])
         Rmse_θΨ_Aver      = Statistics.mean(hydroOther.Rmse_θΨ[1:NiZ])
         Rmse_KΨ_Aver      = Statistics.mean(hydroOther.Rmse_KΨ[1:NiZ])
				
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


end  # module hydrolabOpt
# ............................................................