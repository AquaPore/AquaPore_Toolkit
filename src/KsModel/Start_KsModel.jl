# =============================================================
#		module: startKsModel
# =============================================================

include("θψ_2_KsModel.jl")
include("Opt_KsModel.jl")

module startKsModel
	import ..kunsat, ..optKsModel, ..plot, ..stats, ..θψ_2_KsψModel
	import Statistics
	export START_KSΨMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_KSΨMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_KSΨMODEL(hydro, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param, path; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], Ks_Impermeable=[])
			# NUMBER OF CLASSES
				ClassBool, ClassBool_All, N_Class = KSΨMODEL_CLASS(hydro, NiZ, option, param)

			# DERIVING OBSERVED K(Ψ₁₀ₖₚₐ)
				KΨ_Obs₁₀ₖₚₐ = fill(0.0::Float64, NiZ)
				KΨ_Sim₁₀ₖₚₐ = fill(0.0::Float64, NiZ)
				for iZ=1:NiZ
					KΨ_Obs₁₀ₖₚₐ[iZ] = kunsat.Ψ_2_KUNSAT(option.hydro, 1000.0, iZ, hydro)
				end

			# PERFORM OPTIMISATION OF KsModel ====
			if sum(optimKsmodel.NparamOpt) ≥ 1 && option.data.Kθ && "Ks" ∈ optim.ParamOpt # For security
				for ipClass=1:N_Class
					ClassBool_Select = ClassBool[1:NiZ, ipClass]

					if sum(ClassBool[1:NiZ, ipClass]) ≥ 1 && optimKsmodel.NparamOpt[ipClass] ≥ 1
						println("\n       === ipClass=$ipClass === \n")
						
						KₛModel = optKsModel.START_OPT_KθMODEL(ClassBool_Select, hydro, ipClass, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param)

						# Computing KΨ_Sim₁₀ₖₚₐ
						for iZ=1:NiZ
							if ClassBool_Select[iZ]
								KΨ_Sim₁₀ₖₚₐ[iZ] = θψ_2_KsψModel.KSΨMODEL_START(hydro, ipClass, iZ, ksmodelτ, option, 1000.0; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])
							end
						end

						ksmodelτ = STATISTICS_KSMODEL(ClassBool_Select, hydro, ipClass, KₛModel, ksmodelτ, KΨ_Obs₁₀ₖₚₐ, KΨ_Sim₁₀ₖₚₐ, optimKsmodel, option)

						# All data
						if option.ksModel.Class && option.ksModel.Plot_KsModel
							NameSim = "Class σ_" * string(ipClass)			
							plot.ksmodel.KSMODEL(KₛModel[ClassBool_Select], KΨ_Obs₁₀ₖₚₐ[ClassBool_Select], KΨ_Sim₁₀ₖₚₐ[ClassBool_Select], hydro.Ks[ClassBool_Select], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[ClassBool_Select], hydro.θsMacMat[ClassBool_Select], hydro.σ[ClassBool_Select], option)
						end # if option.ksModel.Class && option.ksModel.Plot_KsModel
					else
						println("\n       === Skipping ipClass=$ipClass === \n")
					end # if optimKsmodel.NparamOpt[ipClass] ≥ 1
				end # for ipClass=1:N_Class

				# Final statistic of all combined Ks data
				if option.ksModel.OptIndivSoil || option.ksModel.Class					
					ksmodelτ = STATISTICS_KSMODEL(ClassBool_All, hydro, 0, KₛModel, ksmodelτ, KΨ_Obs₁₀ₖₚₐ, KΨ_Sim₁₀ₖₚₐ, optimKsmodel, option)
				end

				# PLOTTING ALL SOILS
				if option.ksModel.Plot_KsModel
					NameSim = "All soils"
					plot.ksmodel.KSMODEL(KₛModel[1:NiZ], KΨ_Obs₁₀ₖₚₐ[1:NiZ], KΨ_Sim₁₀ₖₚₐ[1:NiZ], hydro.Ks[1:NiZ], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[1:NiZ], hydro.θsMacMat[1:NiZ], hydro.σ[1:NiZ], option)	
				end

			# RUN KₛModel
			else # ``````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````
				for ipClass=1:N_Class
					ClassBool_Select = ClassBool[1:NiZ, ipClass]

					for iZ=1:NiZ
						if ClassBool_Select[iZ]
							KₛModel[iZ] = θψ_2_KsψModel.KSΨMODEL_START(hydro, ipClass, iZ, ksmodelτ, option, 0.0; Flag_IsTopsoil=false,Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])
							
							KΨ_Sim₁₀ₖₚₐ[iZ] = θψ_2_KsψModel.KSΨMODEL_START(hydro, ipClass, iZ, ksmodelτ, option, 10_00.0; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])
						end
					end # for ipClass=1:N_Class, 
				end #ifor ipClass=1:N_Class

				#Statistics
					ksmodelτ = STATISTICS_KSMODEL(ClassBool_All, hydro, 0, KₛModel, ksmodelτ, KΨ_Obs₁₀ₖₚₐ, KΨ_Sim₁₀ₖₚₐ, optimKsmodel, option)

				# PLOTTING ALL SOILS
				if option.ksModel.Plot_KsModel
					NameSim = "All soils"
					plot.ksmodel.KSMODEL(KₛModel[1:NiZ], KΨ_Obs₁₀ₖₚₐ[1:NiZ], KΨ_Sim₁₀ₖₚₐ[1:NiZ], hydro.Ks[1:NiZ], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[1:NiZ], hydro.θsMacMat[1:NiZ], hydro.σ[1:NiZ], option)	
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
				N_Class = length(param.ksModel.σₛₚₗᵢₜ) - 1
				ClassBool = fill(false::Bool, NiZ, N_Class)

				for ipClass=1:N_Class, iZ=1:NiZ
					σ_Min = param.ksModel.σₛₚₗᵢₜ[ipClass]
					σ_Max = param.ksModel.σₛₚₗᵢₜ[ipClass+1] 

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
		function STATISTICS_KSMODEL(ClassBool, hydro, ipClass, KₛModel, ksmodelτ, KΨ_Obs₁₀ₖₚₐ, KΨ_Sim₁₀ₖₚₐ, optimKsmodel, option)
	
			# STATISTICS
				# For observed and simulated Ks
               Nse_τ₀    = stats.NSE(log10.(hydro.Ks[ClassBool]) , log10.(KₛModel[ClassBool]))
               Rmse_τ₀   = stats.RMSE(log10.(hydro.Ks[ClassBool]) , log10.(KₛModel[ClassBool]))
               σ_τ₀      = Statistics.std(log10.(hydro.Ks[ClassBool]) .- log10.(KₛModel[ClassBool]))
               Wilmot_τ₀ = stats.NSE_WILMOT(log.(hydro.Ks[ClassBool]) , log.(KₛModel[ClassBool]))
               Ccc_τ₀    = stats.stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(log.(hydro.Ks[ClassBool]), log.(KₛModel[ClassBool]))

				# For observed and simulated K(Ψ₁₀ₖₚₐ)
               Nse_KΨ₁₀ₖₚₐ    = stats.NSE(log.(KΨ_Obs₁₀ₖₚₐ[ClassBool]) , log.(KΨ_Sim₁₀ₖₚₐ[ClassBool]))
               Rmse_KΨ₁₀ₖₚₐ   = stats.RMSE(log10.(KΨ_Obs₁₀ₖₚₐ[ClassBool]) , log10.(KΨ_Sim₁₀ₖₚₐ[ClassBool]))
               σ_KΨ₁₀ₖₚₐ      = Statistics.std(log10.(KΨ_Sim₁₀ₖₚₐ[ClassBool]).-log10.(KΨ_Obs₁₀ₖₚₐ[ClassBool]))
               Wilmot_KΨ₁₀ₖₚₐ = stats.NSE_WILMOT(log.(KΨ_Obs₁₀ₖₚₐ[ClassBool]) , log.(KΨ_Sim₁₀ₖₚₐ[ClassBool]))
               Ccc_KΨ₁₀ₖₚₐ    = stats.stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(log.(KΨ_Obs₁₀ₖₚₐ[ClassBool]) , log.(KΨ_Sim₁₀ₖₚₐ[ClassBool]))

				if ipClass == 0
					println("\n       === Statistics all data === \n")
				end

			# PRINING RESULTS
				if hydro.Ks[1] > 0.0
               println("		 Nse_τ          = $(Nse_τ₀)")
               println("		 Rmse_τ         = $(Rmse_τ₀)")
               println("		 σ_τ₀           = $(σ_τ₀)")
               println("		 Wilmot_τ       = $(Wilmot_τ₀)")
               println("		 Ccc_τ          = $(Ccc_τ₀) \n")

               println("		 Nse_KΨ₁₀ₖₚₐ    = $(Nse_KΨ₁₀ₖₚₐ)")
               println("		 Rmse_KΨ₁₀ₖₚₐ   = $(Rmse_KΨ₁₀ₖₚₐ)")
               println("		 σ_KΨ₁₀ₖₚₐ      = $(σ_KΨ₁₀ₖₚₐ)")
               println("		 Wilmot_KΨ₁₀ₖₚₐ = $(Wilmot_KΨ₁₀ₖₚₐ)")
               println("		 Ccc_KΨ₁₀ₖₚₐ    = $(Ccc_KΨ₁₀ₖₚₐ) \n")
				end
				
			if ipClass > 0
				for iParam = 1:optimKsmodel.NparamOpt[ipClass]	
					# Getting the current values of every layer of the hydro parameter of interest
						vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[ipClass, iParam]))
						println("		", Symbol(optimKsmodel.ParamOpt[ipClass, iParam]) , "=" ,vectParam[ipClass])
				end # for loop

				ksmodelτ.Nse_τ[ipClass]    = Nse_τ₀
				ksmodelτ.Rmse_τ[ipClass]   = Rmse_τ₀
				ksmodelτ.Wilmot_τ[ipClass] = Wilmot_τ₀
				ksmodelτ.Ccc_τ[ipClass]    = Ccc_τ₀
			end	
		return ksmodelτ
		end  # function: STATISTICS_KSMODEL
	# ------------------------------------------------------------------

end  # module: startKsModel
# =====================================================================