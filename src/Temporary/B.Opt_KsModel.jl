# =============================================================
#		module: optKsModel
# =============================================================
module optKsModel
	import ..cst, ..kunsat, ..stats, ..θψ_2_KsψModel
	import BlackBoxOptim
	export START_OPT_KθMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_OPT_KθMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_OPT_KθMODEL(ClassBool_Select, hydro, ipClass, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param)
				
			# Deriving the feasible range of the τ parameters
				SearchRange = SEARCHRANGE(ipClass, optimKsmodel)

			# Optimisation algorithme, MaxFuncEvals=1000
				Optimization = BlackBoxOptim.bboptimize(X -> OF_KθMODEL(ClassBool_Select, hydro, ipClass, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[]); SearchRange=SearchRange, NumDimensions=optimKsmodel.NparamOpt[ipClass], TraceMode=:silent)

			# Deriving the optimal τ parameters from X
				X = BlackBoxOptim.best_candidate(Optimization)

			# Putting X parameters into τ
				ksmodelτ = X_2_τ(ipClass, ksmodelτ, optimKsmodel, X)

			# Computing optimal KₛModel
				KₛModel = θψ_2_KsψModel.KSMODEL(ClassBool_Select, hydro, ipClass, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], ipClass=ipClass)

		return KₛModel
		end  # function: START_OPT_KθMODEL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_KθMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_KθMODEL(ClassBool_Select::Vector{Bool}, hydro, ipClass, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], KsMinMax=0.005555556)

			# Deriving the optimal τ parameters from X
				ksmodelτ = X_2_τ(ipClass, ksmodelτ, optimKsmodel, X)

			#	Compuring Ks model
				KₛModel = θψ_2_KsψModel.KSMODEL(ClassBool_Select, hydro, ipClass, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

			if option.ksModel.Opt_Kθ
				Ψ_Obs =  [0.0, 10.0, 20.0, 50.0, 100.0, 500.0, 1000.0, 2000.0, 3300.0, 4000.0, 5000.0,100_00.0, 500_00.0, 1000_00.0]::Float64 # mm
			else
				Ψ_Obs =  [0.0]::Float64 # mm
			end
			N_ΨObs = length(Ψ_Obs)


			Kθ_Log_Obs = fill(0.0::Float64, N_ΨObs)
			Kθ_Log_Sim = fill(0.0::Float64, N_ΨObs)

			for iZ=1:NiZ
				for iΨ =1:N_ΨObs
	 				Kθ_Sim = θψ_2_KsψModel.KSMODEL(ClassBool_Select, hydro, ipClass, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

					Kθ_Log_Sim[iΨ] = log1p(cst.MmS_2_MmH * Kθ_Log_Obs)

					Kθ_Obs = kunsat.Ψ_2_KUNSAT(option.hydro, Ψ_Obs[iΨ], iZ, hydro)

					Kθ_Log_Obs[iΨ] = log1p(cst.MmS_2_MmH * Kθ_Obs)

					Of_Ks = 1.0 - stats.NSE_WILMOT(Kθ_Log_Sim, Kθ_Log_Obs)

				end # iΨ
			end # iZ

			#	Compuring Ks model
				KₛModel = θψ_2_KsψModel.KSMODEL(ClassBool_Select, hydro, ipClass, KₛModel, KₛModel⍰, ksmodelτ, NiZ, optim, optimKsmodel, option, param; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

			# Computing the Objective function ========================================================
				Ks_ObsTransformed = fill(0.0::Float64, NiZ)
				Ks_SimTransformed = fill(0.0::Float64, NiZ)
			
			# Determening when Ks < KsMinMax
				if option.ksModel.Of_Split_KsSlowKsFast
					KsSmall_True = fill(false, NiZ)
					KsLarge_True = fill(false, NiZ)

					for iZ=1:NiZ
						if ClassBool_Select[iZ] # If we have selected the data
							if hydro.Ks[iZ] ≥ KsMinMax
								KsSmall_True[iZ] = false
								KsLarge_True[iZ] = true
							else
								KsSmall_True[iZ] = true
								KsLarge_True[iZ] = false
							end # if hydro.Ks[iZ] ≥ KsMinMax
						else
							KsSmall_True[iZ] = false
							KsLarge_True[iZ] = false
						end # if ClassBool_Select[iZ]
					end # for iZ=1:NiZ
				end # if option.ksModel.Of_Split_KsSlowKsFast
				#____________________________________

			Ks_ObsTransformed = fill(0.0::Float64, NiZ)
			Ks_SimTransformed = fill(0.0::Float64, NiZ)

			# CONVERT UNITS mms-> mm/h
				for iZ=1:NiZ
					Ks_ObsTransformed[iZ] = cst.MmS_2_MmH .* hydro.Ks[iZ]
					Ks_SimTransformed[iZ] = cst.MmS_2_MmH .* KₛModel[iZ]
				end # for iZ=1:NiZ

			if option.ksModel.Of_Split_KsSlowKsFast
				Ks_ObsTransformed_Small = Ks_ObsTransformed[KsSmall_True[1:NiZ]]
				Ks_SimTransformed_Small = Ks_SimTransformed[KsSmall_True[1:NiZ]]

				Ks_ObsTransformed_Large = Ks_ObsTransformed[KsLarge_True[1:NiZ]]
				Ks_SimTransformed_Large = Ks_SimTransformed[KsLarge_True[1:NiZ]]
			end

			# LOG TRANSFORMATION
				if option.ksModel.Of_Split_KsSlowKsFast
					Ks_ObsTransformed_Large = log1p.(Ks_ObsTransformed_Large)
					Ks_SimTransformed_Large = log1p.(Ks_SimTransformed_Large)

				else
					Ks_ObsTransformed = log1p.(Ks_ObsTransformed)
					Ks_SimTransformed = log1p.(Ks_SimTransformed)
				end

				if option.ksModel.Of_Split_KsSlowKsFast
					Of_KsSmall = 1.0 - stats.NSE_WILMOT(Ks_ObsTransformed_Small, Ks_SimTransformed_Small)
					Of_KsLarge = 1.0 - stats.NSE_WILMOT(Ks_ObsTransformed_Large , Ks_SimTransformed_Large)

					Of_Ks = param.ksModel.WeightKsSlow * Of_KsSmall + (1 - param.ksModel.WeightKsSlow) * Of_KsLarge

				else
					Of_Ks = 1.0 - stats.NSE_WILMOT(Ks_ObsTransformed[ClassBool_Select] , Ks_SimTransformed[ClassBool_Select])

				end
	
		return Of_Ks
		end  # function: OF_KSMODELa
		# --------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SEARCHRANGE
	#		Required by BlackBoxOptim
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SEARCHRANGE(ipClass, optimKsmodel)
			ParamOpt_Min₂ = copy(optimKsmodel.ParamOpt_Min[ipClass, 1:optimKsmodel.NparamOpt[ipClass]])
			ParamOpt_Max₂ = copy(optimKsmodel.ParamOpt_Max[ipClass, 1:optimKsmodel.NparamOpt[ipClass]])
		return SearchRange = (collect(zip(Float64.(ParamOpt_Min₂), Float64.(ParamOpt_Max₂))))
		end  # function: SEARCHRANGE
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function X_2_τ(ipClass, ksmodelτ, optimKsmodel, X)
			for iParam = 1:optimKsmodel.NparamOpt[ipClass]
				Paramₐ = X[iParam]
				
				# Getting the current values of every layer of the hydro parameter of interest
					vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[ipClass, iParam]))

				# Updating the value of the parameters for the layer wanting to optimize by keeping the other values constant
					vectParam[ipClass] = Paramₐ

				# Putting the updated hydro into ksmodelτ
					setfield!(ksmodelτ, Symbol(optimKsmodel.ParamOpt[ipClass, iParam]), vectParam)
			end # for loop
		return ksmodelτ
		end  # function: PARAM
	#..................................................................


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : name
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_Kθ(option, hydro, NiZ )


			
		return
		end  # function: name
# ------------------------------------------------------------------
	
end  # module: optKsModel
# ========================================================================