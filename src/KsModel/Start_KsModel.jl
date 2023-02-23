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
			ClassBool, ClassBool_All, N_Class = KSΨMODEL_CLASS(hydro, NiZ, option, param)

			# PERFORM OPTIMISATION OF KsModel ====
			if sum(optimKsmodel.NparamOpt) ≥ 1 && option.data.Kθ && "Ks" ∈ optim.ParamOpt # For security
				for ipClass=1:N_Class

					if optimKsmodel.NparamOpt[ipClass] ≥ 1

						println("\n       === ipClass=$ipClass === \n")
						ClassBool_Select = ClassBool[1:NiZ, ipClass]

						KₛModel = optKsModel.START_OPT_KθMODEL(ClassBool_Select, hydro, ipClass, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param)

						if !(option.ksModel.OptIndivSoil)
							ksmodelτ = STATISTICS_KSMODEL(ClassBool_Select, hydro, ipClass, KₛModel, ksmodelτ, optimKsmodel, option)
						end

						# All data
						if option.ksModel.Class && option.ksModel.Plot_KsModel
							NameSim = "Class σ_" * string(ipClass)
				
							plot.ksmodel.KSMODEL(KₛModel[ClassBool_Select], hydro.Ks[ClassBool_Select], NameSim,path.plotSoilwater.Plot_KsModel, hydro.θr[ClassBool_Select], hydro.θsMacMat[ClassBool_Select], hydro.σ[ClassBool_Select], option)
						end # if option.ksModel.Class && option.ksModel.Plot_KsModel
					end # if optimKsmodel.NparamOpt[ipClass] ≥ 1
				end # for ipClass=1:N_Class


				if option.ksModel.Class
					ksmodelτ = STATISTICS_KSMODEL(ClassBool_All, hydro, 0, KₛModel, ksmodelτ, optimKsmodel, option; 🎏allClass=true)
				end

				# PLOTTING ALL SOILS
				if option.ksModel.Plot_KsModel
					NameSim = "All soils"
					plot.ksmodel.KSMODEL(KₛModel[1:NiZ], hydro.Ks[1:NiZ], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[1:NiZ], hydro.θsMacMat[1:NiZ], hydro.σ[1:NiZ], option)	
				end

			# RUN KₛModel
			else
				for ipClass=1:N_Class
					ClassBool_Select = ClassBool[1:NiZ, ipClass]

					for iZ=1:NiZ
						if ClassBool_Select[iZ]
							KₛModel[iZ] = θψ_2_KsψModel.KSΨMODEL_START(hydro, ipClass, iZ, ksmodelτ, option, 0.0; Flag_IsTopsoil=false,Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])
						end
					end # for ipClass=1:N_Class, 
				end #iZ=1:NiZ

				#Statistics
					ksmodelτ = STATISTICS_KSMODEL(ClassBool_All, hydro, 0, KₛModel, ksmodelτ, optimKsmodel, option; 🎏allClass=true)

				# PLOTTING ALL SOILS
				if option.ksModel.Plot_KsModel
					NameSim = "All soils"
					plot.ksmodel.KSMODEL(KₛModel[1:NiZ], hydro.Ks[1:NiZ], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[1:NiZ], hydro.θsMacMat[1:NiZ], hydro.σ[1:NiZ], option)	
				end
			end  # if: optimKsmodel

			for iZ=1:NiZ
				if "Ks" ∉ optim.ParamOpt
					hydro.Ks[iZ] = KₛModel[iZ]

					# If wanting to assure that the feasible range is physical
					hydro.Ks[iZ] = max( min(hydro.Ks[iZ], hydro.Ks_Max[iZ]), hydro.Ks_Min[iZ])
				end #  hydro.Ks[iZ] < eps(100.0)
			end # if: hydro.Ks[iZ] > eps(10.0)


			# SMAP SPECIAL COORECTING FOR IMPERMEABLE LAYERS
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

			ClassBool_All = fill(true::Bool, NiZ)

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
				
		return ClassBool, ClassBool_All, N_Class
		end  # function: SELECTION
	# ------------------------------------------------------------------
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : STATISTICS_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function STATISTICS_KSMODEL(ClassBool_Select, hydro, ipClass, KₛModel, ksmodelτ, optimKsmodel, option; 🎏allClass=false)
			# STATISTICS
            Nse_τ₀    = stats.NSE(log1p.(hydro.Ks[ClassBool_Select]) , log1p.(KₛModel[ClassBool_Select]))
            Rmse_τ₀   = stats.RMSE(log1p.(hydro.Ks[ClassBool_Select]) , log1p.(KₛModel[ClassBool_Select]))
            Wilmot_τ₀ = stats.NSE_WILMOT(log1p.(hydro.Ks[ClassBool_Select]) , log1p.(KₛModel[ClassBool_Select]))
            Ccc_τ₀    = stats.stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(log1p.(hydro.Ks[ClassBool_Select]) , log1p.(KₛModel[ClassBool_Select]))

				if !(🎏allClass)
               ksmodelτ.Nse_τ[ipClass]    = Nse_τ₀
               ksmodelτ.Rmse_τ[ipClass]   = Rmse_τ₀
               ksmodelτ.Wilmot_τ[ipClass] = Wilmot_τ₀
               ksmodelτ.Ccc_τ[ipClass]    = Ccc_τ₀
				end  # if: !(🎏allClass)

				if 🎏allClass
					println("\n       === Statistics all data === \n")
				end
			# PRINING RESULTS
				if hydro.Ks[1] > 1
					println("		 Nse_τ    =  $(Nse_τ₀)")
					println("		 Rmse_τ   =  $(Rmse_τ₀)")
					println("		 Wilmot_τ =  $(Wilmot_τ₀)")
					println("		 Ccc_τ    =  $(Ccc_τ₀) \n")
				end

				if !(🎏allClass)
				for iParam = 1:optimKsmodel.NparamOpt[ipClass]	
					# Getting the current values of every layer of the hydro parameter of interest
						vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[ipClass, iParam]))

						println("		", Symbol(optimKsmodel.ParamOpt[ipClass, iParam]) , "=" ,vectParam)
				end # for loop
				end # 🎏allClass
			
		return ksmodelτ
		end  # function: STATISTICS_KSMODEL
	# ------------------------------------------------------------------

end  # module: startKsModel
# =====================================================================