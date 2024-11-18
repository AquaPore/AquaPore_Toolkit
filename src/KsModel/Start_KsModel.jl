# =============================================================
#		module: startKsModel
# =============================================================

include("θψ_2_KsModel.jl")
include("Opt_KsModel.jl")

module startKsModel
	import ..kunsat, ..optKsModel, ..plot, ..stats, ..θψ_2_KsψModel, ..cst
	import Statistics
	export START_KSΨMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_KSΨMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_KSΨMODEL(hydro, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param, path; IsTopsoil=[], RockFragment=[], Ks_Impermeable=[], ∑Psd=[])
			# Plotting options
				🎏_Plot_Rf     = false
				🎏_Plot_KsFunc = false
				🎏_Plot_KsFunc = false
				🎏_Plot_Tclay  = false
				
			# WHAT DATA DO WE HAVE FOR OUR ANALYSIS 
				🎏_Clay         = !isempty(∑Psd)
				🎏_IsTopsoil    = !isempty(IsTopsoil)
				🎏_RockFragment = !isempty(RockFragment) && option.run.RockCorection

			# TIME NOW 
				Time_Start = time()

			# SEPERATING DATA INTO CLASSES
				ClassBool, ClassBool_All, N_Class = startKsModel.KSΨMODEL_CLASS(NiZ, optimKsmodel, option)

			# DERIVING OBSERVED K(Ψ₁₀ₖₚₐ)
				KΨ_Obs₁₀ₖₚₐ = fill(0.0::Float64, NiZ)
				KΨ_Sim₁₀ₖₚₐ = fill(0.0::Float64, NiZ)
				for iZ=1:NiZ
					KΨ_Obs₁₀ₖₚₐ[iZ] = kunsat.KUNSAT_θΨSe(option.hydro, 10_00.0, iZ, hydro)
				end

			# OPTIMISATION OR RUNNING
				for ipClass=1:N_Class
					# Selecting the data which contains the class of interest
					ClassBool_Select = ClassBool[1:NiZ, ipClass]
					
					# Do we need to optimise the class
					if optimKsmodel.NparamOpt[ipClass] ≥ 1 && option.data.Kθ && "Ks" ∈ optim.ParamOpt # For security to determine if we have 
						printstyled("\n       === ipClass=$ipClass === \n"; color=:yellow)

						# Optimising the model	
							KₛModel = optKsModel.START_OPT_KθMODEL(∑Psd, 🎏_Clay, ClassBool_Select, hydro, ipClass, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param; 🎏_IsTopsoil=🎏_IsTopsoil, 🎏_RockFragment=🎏_RockFragment, RockFragment=RockFragment, IsTopsoil=IsTopsoil)

						# Plotting
							if option.ksModel.Plot_KsModel
								NameSim = "Class_" * string(ipClass)			
								plot.ksmodel.KSMODEL(KₛModel[ClassBool_Select], KΨ_Obs₁₀ₖₚₐ[ClassBool_Select], KΨ_Sim₁₀ₖₚₐ[ClassBool_Select], hydro.Ks[ClassBool_Select], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[ClassBool_Select], hydro.θsMacMat[ClassBool_Select], hydro.σ[ClassBool_Select], option)
							end # option.ksModel.Plot_KsModel
					# ~~~~~~~~~~~~~~~~
					else
					# ~~~~~~~~~~~~~~~~ 			
						for iZ=1:NiZ
							if ClassBool_Select[iZ]
								KₛModel[iZ] = θψ_2_KsψModel.KSΨMODEL_START(∑Psd, 🎏_Clay, hydro, ipClass, iZ, ksmodelτ, option, param, 0.0; 🎏_IsTopsoil=🎏_IsTopsoil,🎏_RockFragment=🎏_RockFragment, IsTopsoil=IsTopsoil, RockFragment=RockFragment)
							end # if ClassBool_Select[iZ]
						end # for iZ=1:NiZ 

					end # optimKsmodel.NparamOpt[ipClass] ≥ 1	

					# Computing KΨ_Sim₁₀ₖₚₐ
						for iZ=1:NiZ
							if ClassBool_Select[iZ]
								KΨ_Sim₁₀ₖₚₐ[iZ] = θψ_2_KsψModel.KSΨMODEL_START(∑Psd, 🎏_Clay, hydro, ipClass, iZ, ksmodelτ, option, param, 1000.0; 🎏_IsTopsoil=🎏_IsTopsoil, 🎏_RockFragment=🎏_RockFragment, RockFragment=RockFragment, IsTopsoil=IsTopsoil)
							end
						end
					
					if option.data.Kθ
						ksmodelτ = startKsModel.STATISTICS_KSMODEL(ClassBool_Select, hydro, ipClass, KₛModel, ksmodelτ, KΨ_Obs₁₀ₖₚₐ, KΨ_Sim₁₀ₖₚₐ, optimKsmodel, option)
					end	
				end # for ipClass=1:N_Class

				# ======================================================================================================
				# 				FOR ALL SOILS
				# ======================================================================================================

					printstyled("\n       === ~ALL SOILS~ === \n", color=:yellow)

				# CHECKING FOR CONSISTENCY & BOUNDARIES
					for iZ=1:NiZ
						KₛModel[iZ] = max(min(KₛModel[iZ], hydro.Ks_Max[iZ]), hydro.Ks_Min[iZ])
	
						if option.run.Smap
							# Special cases for impermeable layers
							if Ks_Impermeable[iZ] ≥ 0.0
								KₛModel[iZ] = Ks_Impermeable[iZ]
							end
						end
	
						if "Ks" ∉ optim.ParamOpt
							hydro.Ks[iZ] = KₛModel[iZ]
						end #  hydro.Ks[iZ] < eps(100.0)
					end
			
				# # STATISTICS
					if option.data.Kθ				
				 		~ = STATISTICS_KSMODEL(ClassBool_All, hydro, 0, KₛModel, ksmodelτ, KΨ_Obs₁₀ₖₚₐ, KΨ_Sim₁₀ₖₚₐ, optimKsmodel, option)
					end

			# PLOTTING ALL SOILS
				if option.ksModel.Plot_KsModel && option.data.Kθ
					NameSim = "Allsoils"
					plot.ksmodel.KSMODEL(KₛModel[1:NiZ], KΨ_Obs₁₀ₖₚₐ[1:NiZ], KΨ_Sim₁₀ₖₚₐ[1:NiZ], hydro.Ks[1:NiZ], NameSim, path.plotSoilwater.Plot_KsModel, hydro.θr[1:NiZ], hydro.θsMacMat[1:NiZ], hydro.σ[1:NiZ], option)

					if 🎏_Plot_Tclay
						plot.ksmodel.KSMODEL_TCLAY( path.plotSoilwater.Plot_KsModel, option, ksmodelτ, 1)
					end
					if 🎏_Plot_KsFunc
						plot.ksmodel.KSMODEL_FUNCTIONS( path.plotSoilwater.Plot_KsModel, option, ksmodelτ, 1)
					end
					if 🎏_Plot_Rf
						plot.ksmodel.KSMODEL_RF( path.plotSoilwater.Plot_KsModel, hydro, option, ksmodelτ, 1)
					end
				end

			Time_End = time()
			println("\n		~~~~~~~~~~ Time of simulations $(floor(Time_End-Time_Start)) Seconds \n")
		return hydro, KₛModel, N_Class
		end  # function: START_KSΨMODEL
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSΨMODEL_CLASS 
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSΨMODEL_CLASS(NiZ, optimKsmodel, option)
			ClassBool_All = fill(true::Bool, NiZ)

			if option.ksModel.OptIndivSoil
				N_Class = NiZ 
				ClassBool = fill(false::Bool, NiZ, N_Class)
				for iZ=1:NiZ
					ClassBool[iZ, iZ] = true
				end

			else
				N_Class = optimKsmodel.N_KsClass
				ClassBool = fill(false::Bool, NiZ, N_Class)

				for ipClass=1:N_Class, iZ=1:NiZ
					if optimKsmodel.KsClass[iZ] == ipClass
						ClassBool[iZ, ipClass] = true
					end
				end
			end
				
		return ClassBool, ClassBool_All, N_Class
		end  # function: SELECTION
	# ------------------------------------------------------------------
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : STATISTICS_KSMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function STATISTICS_KSMODEL(ClassBool_Select, hydro, ipClass, KₛModel, ksmodelτ, KΨ_Obs₁₀ₖₚₐ, KΨ_Sim₁₀ₖₚₐ, optimKsmodel, option)
	
			# For observed and simulated Ks
            Nse_τ₀    = stats.NSE(log10.(cst.MmS_2_MmH * hydro.Ks[ClassBool_Select]) , log10.(cst.MmS_2_MmH * KₛModel[ClassBool_Select]))
            Rmse_τ₀   = stats.RMSE( log10.(cst.MmS_2_MmH * hydro.Ks[ClassBool_Select]) , log10.(cst.MmS_2_MmH * KₛModel[ClassBool_Select]))
            σ_τ₀      = Statistics.std(log10.(cst.MmS_2_MmH * hydro.Ks[ClassBool_Select]) .- log10.(cst.MmS_2_MmH * KₛModel[ClassBool_Select]))
            Wilmot_τ₀ = stats.NSE_WILMOT(log10.(cst.MmS_2_MmH * hydro.Ks[ClassBool_Select]) , log10.(cst.MmS_2_MmH * KₛModel[ClassBool_Select]))
            Ccc_τ₀    = stats.stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(log10.(cst.MmS_2_MmH .* hydro.Ks[ClassBool_Select]), log10.(cst.MmS_2_MmH .* KₛModel[ClassBool_Select]))

			# For observed and simulated K(Ψ₁₀ₖₚₐ)
            Nse_KΨ₁₀ₖₚₐ    = stats.NSE(log10.(cst.MmS_2_MmH * KΨ_Obs₁₀ₖₚₐ[ClassBool_Select]) , log10.(cst.MmS_2_MmH * KΨ_Sim₁₀ₖₚₐ[ClassBool_Select]))
            Rmse_KΨ₁₀ₖₚₐ   = stats.RMSE(log10.(cst.MmS_2_MmH * KΨ_Obs₁₀ₖₚₐ[ClassBool_Select]) , log10.(cst.MmS_2_MmH * KΨ_Sim₁₀ₖₚₐ[ClassBool_Select]))
            σ_KΨ₁₀ₖₚₐ      = Statistics.std(log10.(KΨ_Sim₁₀ₖₚₐ[ClassBool_Select]).-log10.(KΨ_Obs₁₀ₖₚₐ[ClassBool_Select]))
            Wilmot_KΨ₁₀ₖₚₐ = stats.NSE_WILMOT(log10.(cst.MmS_2_MmH * KΨ_Obs₁₀ₖₚₐ[ClassBool_Select]) , log10.(cst.MmS_2_MmH * KΨ_Sim₁₀ₖₚₐ[ClassBool_Select]))
            Ccc_KΨ₁₀ₖₚₐ    = stats.stats.NSE_CONCORDANCE_CORELATION_COEFICIENT(log10.(cst.MmS_2_MmH * KΨ_Obs₁₀ₖₚₐ[ClassBool_Select]) , log10.(cst.MmS_2_MmH * KΨ_Sim₁₀ₖₚₐ[ClassBool_Select]))

			# PRINING RESULTS
				println("		 Nse_τ          = $(Nse_τ₀) log10 mm/h")
				println("		 Rmse_τ         = $(Rmse_τ₀) log10 mm/h")
				println("		 σ_τ₀           = $(σ_τ₀) log10 mm/h")
				println("		 Wilmot_τ       = $(Wilmot_τ₀) log10 mm/h")
				println("		 Ccc_τ          = $(Ccc_τ₀) log10 mm/h \n")

				println("		 Nse_KΨ₁₀ₖₚₐ    = $(Nse_KΨ₁₀ₖₚₐ)")
				println("		 Rmse_KΨ₁₀ₖₚₐ   = $(Rmse_KΨ₁₀ₖₚₐ)")
				println("		 σ_KΨ₁₀ₖₚₐ      = $(σ_KΨ₁₀ₖₚₐ)")
				println("		 Wilmot_KΨ₁₀ₖₚₐ = $(Wilmot_KΨ₁₀ₖₚₐ)")
				println("		 Ccc_KΨ₁₀ₖₚₐ    = $(Ccc_KΨ₁₀ₖₚₐ) \n")
				
			if ipClass ≥ 1
				for iParam = 1:optimKsmodel.NparamOpt[ipClass]	
					# Getting the current values of every layer of the hydro parameter of interest
						vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[ipClass, iParam]))
						println("		", Symbol(optimKsmodel.ParamOpt[ipClass, iParam]) , " = " ,vectParam[ipClass])
				end # for loop

				ksmodelτ.Nse_τ[ipClass]    = Nse_τ₀
				ksmodelτ.Rmse_τ[ipClass]   = Rmse_τ₀
				ksmodelτ.Wilmot_τ[ipClass] = Wilmot_τ₀
				ksmodelτ.Ccc_τ[ipClass]    = Ccc_τ₀
			end # if ipClass ≥ 1	
		return ksmodelτ
		end  # function: STATISTICS_KSMODEL
	# ------------------------------------------------------------------

end  # module: startKsModel
# =====================================================================