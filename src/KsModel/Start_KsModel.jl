# =============================================================
#		module: startKsModel
# =============================================================

include("θψ_2_KsModel.jl")
include("Opt_KsModel.jl")

module startKsModel
	import ..θψ2KsModel, ..optKsModel, ..stats, ..plot
	export START_KθMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_KθMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_KθMODEL(hydro, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param, path; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], Ks_Impermeable=[], ipLayer=1, N_Group=2)

			# GROUPING CLASSES
				GroupBool, ipGroup, N_Group = GROUPING_KSMODEL(hydro, N_Group, NiZ, option, param)

			# If there are τ parameters to be optimised
			if sum(optimKsmodel.NparamOpt) ≥ 1 && option.data.Kθ && "Ks" ∈ optim.ParamOpt # For security
				for iGroup_Opt=1:N_Group

					println("\n       === iGroup_Opt=$iGroup_Opt === \n")

					# Quick fix
					if optimKsmodel.NparamOpt[iGroup_Opt] ≥ 1
						GroupBool_Select = GroupBool[1:NiZ, iGroup_Opt]

						KₛModel = optKsModel.START_OPT_KθMODEL(GroupBool_Select, hydro, iGroup_Opt, ipGroup, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param)

						ksmodelτ = STATISTICS_KSMODEL(hydro,iGroup_Opt, GroupBool_Select, KₛModel, ksmodelτ, optimKsmodel, option)

						if option.ksModel.Plot_KsModel && !(option.ksModel.OptIndivSoil)
							NameSim = "σ_" * string(iGroup_Opt)

							plot.ksmodel.KSMODEL(KₛModel[GroupBool_Select], hydro.Ks[GroupBool_Select], NameSim,path.plotSoilwater.Plot_KsModel, hydro.θr[GroupBool_Select], hydro.θsMacMat[GroupBool_Select], hydro.σ[GroupBool_Select], option)
						end # if option.Plot

					end # if optimKsmodel.NparamOpt[iGroup_Opt] ≥ 1
				end # for iGroup_Opt=1:N_Group

				if option.ksModel.Plot_KsModel && option.ksModel.OptIndivSoil
					NameSim = "All_"
					plot.ksmodel.KSMODEL(KₛModel[1:NiZ], hydro.Ks[1:NiZ], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[1:NiZ], hydro.θsMacMat[1:NiZ], hydro.σ[1:NiZ], option)	
				end

			# RUN KₛModel
			else
				GroupBool_Select = fill(true, NiZ)

				iGroup_Opt = 1

				KₛModel = θψ2KsModel.KSMODEL(GroupBool_Select, hydro, ipGroup, KₛModel, ksmodelτ, NiZ::Int64, optim, optimKsmodel, option, param; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

				~ = STATISTICS_KSMODEL(hydro, iGroup_Opt, GroupBool_Select, KₛModel, ksmodelτ, optimKsmodel)
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
		
		return hydro, KₛModel, N_Group
		end  # function: START_KθMODEL
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : GROUPING_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function GROUPING_KSMODEL(hydro, N_Group, NiZ, option, param)

			if option.ksModel.OptIndivSoil
				N_Group = NiZ 
				GroupBool = fill(false::Bool, NiZ, N_Group)
			else
				GroupBool = fill(true::Bool, NiZ, N_Group)
			end
			
         ipGroup   = fill(1::Int64, NiZ)
         
			# Selecting groups if required
			if option.ksModel.Group
				for iZ=1:NiZ

					if !(option.ksModel.OptIndivSoil)
						if hydro.σ[iZ] ≤ param.ksModel.σₛₚₗᵢₜ
							ipGroup[iZ] = 1
							GroupBool[iZ, 1] = true
							GroupBool[iZ, 2] = false
						else
							ipGroup[iZ] = 2
							GroupBool[iZ, 1] = false
							GroupBool[iZ, 2] = true
						end  # if: hydro.

						N_Group = 2
					else
						GroupBool[iZ, iZ] = true
						ipGroup[iZ] = iZ
						N_Group = NiZ
					end  
				end # for iZ=1:NiZ
				
			else
				N_Group = 1
			end #  option.ksModel.Group
				
		return GroupBool, ipGroup, N_Group
		end  # function: SELECTION
	# ------------------------------------------------------------------
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : STATISTICS_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function STATISTICS_KSMODEL(hydro,  iGroup_Opt, GroupBool_Select, KₛModel, ksmodelτ, optimKsmodel, option)
			# STATISTICS
				ksmodelτ.Nse_τ[iGroup_Opt]    = stats.NSE(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))
				ksmodelτ.Rmse_τ[iGroup_Opt]   = stats.RMSE(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))
				ksmodelτ.Wilmot_τ[iGroup_Opt] = stats.NSE_WILMOT(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))
				ksmodelτ.Ccc_τ[iGroup_Opt]    = stats.stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(log1p.(hydro.Ks[GroupBool_Select]) , log1p.(KₛModel[GroupBool_Select]))

			# PRINING RESULTS
				if sum(optimKsmodel.NparamOpt) ≥ 1
					println("		 Nse_τ    =  $(ksmodelτ.Nse_τ)")
					println("		 Rmse_τ   =  $(ksmodelτ.Rmse_τ)")
					println("		 Wilmot_τ =  $(ksmodelτ.Wilmot_τ)")
					println("		 Ccc_τ    =  $(ksmodelτ.Ccc_τ)")
				end

				for iParam = 1:optimKsmodel.NparamOpt[iGroup_Opt]	
					# Getting the current values of every layer of the hydro parameter of interest
						vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[iGroup_Opt, iParam]))

						println("		", Symbol(optimKsmodel.ParamOpt[iGroup_Opt, iParam]) , "=" ,vectParam)
				end # for loop
			
		return ksmodelτ
		end  # function: STATISTICS_KSMODEL
	# ------------------------------------------------------------------


end  # module: startKsModel
# =====================================================================