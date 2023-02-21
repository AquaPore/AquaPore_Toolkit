# =============================================================
#		module: startKsModel
# =============================================================

include("θψ_2_KsModel.jl")
include("Opt_KsModel.jl")

module startKsModel
	import ..optKsModel, ..plot, ..stats, ..θψ_2_KsψModel
	export START_KSΨMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_KSΨMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_KSΨMODEL(hydro, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param, path; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], Ks_Impermeable=[])

			# NUMBER OF CLASSES
				ClassBool, N_Class = KSΨMODEL_CLASS(hydro, NiZ, option, param)

			# If there are τ parameters to be optimised
			if sum(optimKsmodel.NparamOpt) ≥ 1 && option.data.Kθ && "Ks" ∈ optim.ParamOpt # For security
				for ipClass=1:N_Class

					println("\n       === ipClass=$ipClass === \n")

					# Quick fix
					if optimKsmodel.NparamOpt[ipClass] ≥ 1
						GroupBool_Select = ClassBool[1:NiZ, ipClass]

						KₛModel = optKsModel.START_OPT_KθMODEL(GroupBool_Select, hydro, ipClass, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param)

						ksmodelτ = STATISTICS_KSMODEL(hydro,ipClass, GroupBool_Select, KₛModel, ksmodelτ, optimKsmodel, option)

						if option.ksModel.Plot_KsModel && !(option.ksModel.OptIndivSoil)
							NameSim = "σ_" * string(ipClass)

							plot.ksmodel.KSMODEL(KₛModel[GroupBool_Select], hydro.Ks[GroupBool_Select], NameSim,path.plotSoilwater.Plot_KsModel, hydro.θr[GroupBool_Select], hydro.θsMacMat[GroupBool_Select], hydro.σ[GroupBool_Select], option)
						end # if option.Plot

					end # if optimKsmodel.NparamOpt[ipClass] ≥ 1
				end # for ipClass=1:N_Class

				if option.ksModel.Plot_KsModel && option.ksModel.OptIndivSoil
					NameSim = "All_"
					plot.ksmodel.KSMODEL(KₛModel[1:NiZ], hydro.Ks[1:NiZ], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[1:NiZ], hydro.θsMacMat[1:NiZ], hydro.σ[1:NiZ], option)	
				end

			# RUN KₛModel
			else
				GroupBool_Select = fill(true, NiZ)

				ipClass = 1

				KₛModel = θψ_2_KsψModel.KSMODEL(GroupBool_Select, hydro, ipClass, KₛModel, ksmodelτ, NiZ::Int64, optim, optimKsmodel, option, param; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

				~ = STATISTICS_KSMODEL(hydro, ipClass, GroupBool_Select, KₛModel, ksmodelτ, optimKsmodel)
			end  # if: optimKsmodel

			if option.ksModel.Plot_KsModel && option.data.Kθ && "Ks"∈ optim.ParamOpt
				println("\n === ALL SIMULATIONS === \n")
				NameSim = "All_"
				plot.ksmodel.KSMODEL(KₛModel[1:NiZ], hydro.Ks[1:NiZ], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[1:NiZ], hydro.θsMacMat[1:NiZ], hydro.σ[1:NiZ], option)	
			end
			
			for iZ=1:NiZ
				if "Ks" ∉ optim.ParamOpt
					hydro.Ks[iZ] = KₛModel[iZ]

					# If wanting to assure that the feasible range is physical
					hydro.Ks[iZ] = max( min(hydro.Ks[iZ], hydro.Ks_Max[iZ]), hydro.Ks_Min[iZ])
				end #  hydro.Ks[iZ] < eps(100.0)
			end # if: hydro.Ks[iZ] > eps(10.0)


			# COORECTION FOR IMPERMEABLE LAYERS
				if option.run.Smap
					for iZ=1:NiZ
						if Ks_Impermeable[iZ] ≥ 0.0
							KₛModel[iZ] = Ks_Impermeable[iZ]
						end
					end
				end
		
		return hydro, KₛModel, N_Class
		end  # function: START_KSΨMODEL
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSΨMODEL_CLASS 
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSΨMODEL_CLASS(hydro, NiZ, option, param)

			if option.ksModel.OptIndivSoil
				N_Class = NiZ 
				ClassBool = fill(false::Bool, NiZ, N_Class)
				for iZ=1:NiZ
					ClassBool[iZ, iZ] = true
				end

			elseif option.ksModel.Class
				N_Class = length(param.ksModel.σηₛₚₗᵢₜ) - 1
				ClassBool = fill(false::Bool, NiZ, N_Class)

				for ipClass=1:N_Class, iZ=1:NiZ
					σ_Min = param.ksModel.σηₛₚₗᵢₜ[ipClass] * (hydro.σ_Max[iZ] - hydro.σ_Min[iZ]) + hydro.σ_Min[iZ]
					σ_Max = param.ksModel.σηₛₚₗᵢₜ[ipClass+1] * (hydro.σ_Max[iZ] - hydro.σ_Min[iZ]) + hydro.σ_Min[iZ]

					if σ_Min ≤ hydro.σ[iZ] < σ_Max
						ClassBool[iZ, ipClass] = true
					end
				end

			else
				N_Class = 1
				ClassBool = fill(true::Bool, NiZ, N_Class) 
			end
				
		return ClassBool, N_Class
		end  # function: SELECTION
	# ------------------------------------------------------------------
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : STATISTICS_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function STATISTICS_KSMODEL(hydro,  ipClass, GroupBool_Select, KₛModel, ksmodelτ, optimKsmodel, option)
			# STATISTICS
				ksmodelτ.Nse_τ[ipClass]    = stats.NSE(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))
				ksmodelτ.Rmse_τ[ipClass]   = stats.RMSE(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))
				ksmodelτ.Wilmot_τ[ipClass] = stats.NSE_WILMOT(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))
				ksmodelτ.Ccc_τ[ipClass]    = stats.stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))

			# PRINING RESULTS
				if sum(optimKsmodel.NparamOpt) ≥ 1
					println("		 Nse_τ    =  $(ksmodelτ.Nse_τ)")
					println("		 Rmse_τ   =  $(ksmodelτ.Rmse_τ)")
					println("		 Wilmot_τ =  $(ksmodelτ.Wilmot_τ)")
					println("		 Ccc_τ    =  $(ksmodelτ.Ccc_τ)")
				end

				for iParam = 1:optimKsmodel.NparamOpt[ipClass]	
					# Getting the current values of every layer of the hydro parameter of interest
						vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[ipClass, iParam]))

						println("		", Symbol(optimKsmodel.ParamOpt[ipClass, iParam]) , "=" ,vectParam)
				end # for loop
			
		return ksmodelτ
		end  # function: STATISTICS_KSMODEL
	# ------------------------------------------------------------------

end  # module: startKsModel
# =====================================================================