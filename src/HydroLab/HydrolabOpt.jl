# =============================================================
#		module: hypixOpt
# =============================================================
module hydrolabOpt

	import   ..stats, ..optIndivSoil, ..optAllSoil, ..ofHydrolab, ..hydroRelation
	using  Statistics
	export HYDROLABOPT_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIXOPT_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYDROLABOPT_START(;∑Psd, hydro, hydroOther, K_KΨobs=zeros(Float64,1,1), N_KΨobs=1, N_θΨobs, NiZ, optim, optimAllSoils, option, optionₘ, param, θ_θΨobs, Ψ_KΨobs=zeros(Float64,1,1), Ψ_θΨobs)

		# if optionₘ.ΨmacMat_2_σMac_ΨmMac
		# 	for iZ=1:NiZ
		# 		hydro.σMac[iZ]  = hydroRelation.FUNC_ΨmacMat_2_σMac(ΨmacMat=hydro.ΨmacMat[iZ])
					
		# 		hydro.ΨmMac[iZ] = hydroRelation.FUNC_ΨmacMat_2_ΨmMac(ΨmacMat=hydro.ΨmacMat[iZ])
		# 	end
		# end
		

		if optimAllSoils.🎏_Opt
			# OPTIMISATION ALL SOILS
				hydro = optAllSoil.OPTIMIZE_ALLSOILS(;∑Psd, hydro, hydroOther, K_KΨobs, N_KΨobs, N_θΨobs, NiZ, optim, optimAllSoils, option, optionₘ, param, θ_θΨobs, θϵ=0.005, Ψ_KΨobs, Ψ_θΨobs)

				 hydro, hydroOther, Of_Sample = optIndivSoil.OPTIMIZE_INDIVIDUALSOILS(;∑Psd, hydro, hydroOther, K_KΨobs, N_KΨobs, N_θΨobs, NiZ, optim, optimAllSoils, option, optionₘ, param, θ_θΨobs, θϵ=0.005, Ψ_KΨobs, Ψ_θΨobs)

				OfAllSoil = ofHydrolab.OF_ALLSOILS(NiZ::Int64, Of_Sample::Vector{Float64})

				printstyled("				OfAllSoil =  ", trunc(OfAllSoil,digits=3), "\n"; color=:blue)

		else
		# OPTIMISATION INDIVIDUAL SOIL
			hydro, hydroOther, Of_Sample = optIndivSoil.OPTIMIZE_INDIVIDUALSOILS(;∑Psd, hydro, hydroOther, K_KΨobs, N_KΨobs, N_θΨobs, NiZ, optim, optimAllSoils, option, optionₘ, param, θ_θΨobs, θϵ=0.005, Ψ_KΨobs, Ψ_θΨobs)
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
         Nse_θΨ_Aver       = Statistics.mean(max.(hydroOther.Nse_θΨ[1:NiZ],0.0))
         Nse_KΨ_Aver       = Statistics.mean(max.(hydroOther.Nse_KΨ[1:NiZ], 0.0))

         NseWilmot_θΨ_Aver = Statistics.mean(max.(hydroOther.NseWilmot_θΨ[1:NiZ], 0.0))
         NseWilmot_KΨ_Aver = Statistics.mean(max.(hydroOther.NseWilmot_KΨ[1:NiZ], 0.0))

         Rmse_Aver         = Statistics.mean(hydroOther.Rmse[1:NiZ])
         Rmse_θΨ_Aver      = Statistics.mean(hydroOther.Rmse_θΨ[1:NiZ])
         Rmse_KΨ_Aver      = Statistics.mean(hydroOther.Rmse_KΨ[1:NiZ])
				
			if "Ks" ∈ optim.ParamOpt
				Nse_Aver = (Nse_θΨ_Aver + Nse_KΨ_Aver) / 2.0
			else
				Nse_Aver = Nse_θΨ_Aver
			end

			printstyled("	=== === Optimizing Hydraulic parameters  === \n", color=:green)
			printstyled("    		~  Nse_θΨ= $(round(Nse_θΨ_Aver,digits=3)), Nse_KΨ= $(round(Nse_KΨ_Aver,digits=3)), Nse = $(round(Nse_Aver,digits=3))  ~\n", color=:cyan)
			printstyled("    		~  NseWilmot_θΨ= $(round(NseWilmot_θΨ_Aver,digits=3)), NseWilmot_KΨ= $(round(NseWilmot_KΨ_Aver,digits=3)), NseWilmot= $(round((NseWilmot_θΨ_Aver+NseWilmot_KΨ_Aver)*0.5,digits=3)) ~\n", color=:cyan)
			printstyled("    		~  Rmse_θΨ = $(round(Rmse_θΨ_Aver,digits=4)),  RmseLog10_KΨ = $(round(Rmse_KΨ_Aver,digits=4)), Rmse = $(round(Rmse_Aver,digits=4))  ~ \n", color=:cyan)
			# println( "	=== === ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ === ===")

	return hydro, hydroOther
	end  # function: HYPIXOPT_START


end  # module hydrolabOpt
# ............................................................