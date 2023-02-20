# =============================================================
#		module: optKsModel
# =============================================================
module optKsModel
	import ..cst, ..kunsat, ..stats, ..θψ2KsModel
	import BlackBoxOptim
	export START_OPT_KθMODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_OPT_KθMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function START_OPT_KθMODEL(GroupBool_Select, hydro, iGroup_Opt, ipGroup, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param)
				
			# Deriving the feasible range of the τ parameters
				SearchRange = SEARCHRANGE(iGroup_Opt, optimKsmodel)

			# Optimisation algorithme, MaxFuncEvals=1000
				Optimization = BlackBoxOptim.bboptimize(X -> OF_KθMODEL(GroupBool_Select, hydro, iGroup_Opt, ipGroup, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[]); SearchRange=SearchRange, NumDimensions=optimKsmodel.NparamOpt[iGroup_Opt], TraceMode=:silent)

			# Deriving the optimal τ parameters from X
				X = BlackBoxOptim.best_candidate(Optimization)

			# Putting X parameters into τ
				ksmodelτ = X_2_τ(iGroup_Opt, ksmodelτ, optimKsmodel, X)

			# Computing optimal KₛModel
			for iZ=1:NiZ
				KₛModel[iZ] = θψ2KsModel.KθMODEL(hydro, ipGroup, iZ, ksmodelτ, KₛModel, option, 0.0; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], iGroup_Opt=iGroup_Opt)
			end

		return KₛModel
		end  # function: START_OPT_KθMODEL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_KθMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_KθMODEL(GroupBool_Select, hydro, iGroup_Opt, ipGroup, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param, X; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[], KsMinMax=0.005555556)

			# Deriving the optimal τ parameters from X
				ksmodelτ = X_2_τ(iGroup_Opt, ksmodelτ, optimKsmodel, X)

			# If optimising the whole K(Ψ)
				if option.ksModel.Opt_Kθ
					Ψ_Obs =  [0.0, 10.0, 20.0, 50.0, 100.0, 500.0, 1000.0, 2000.0, 3300.0, 4000.0, 5000.0,100_00.0, 500_00.0, 1000_00.0]::Vector{Float64} # mm
				else
					Ψ_Obs =  [0.0]::Vector{Float64} # mm
				end
				N_ΨObs = length(Ψ_Obs)

				Kθ_Log_Obs = fill(0.0::Float64, N_ΨObs)
				Kθ_Log_Sim = fill(0.0::Float64, N_ΨObs)

			# Computing K(Ψ)
			Of_Kθ = 0.0
			for iZ=1:NiZ
			if GroupBool_Select[iZ]
				for iΨ =1:N_ΨObs
					# K(Ψ) simulated
						Kθ_Sim = θψ2KsModel.KθMODEL(hydro, ipGroup, iZ, ksmodelτ, KₛModel, option, Ψ_Obs[iΨ]; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

						Kθ_Log_Sim[iΨ] = log1p(cst.MmS_2_MmH * Kθ_Sim[ipGroup[iZ]])

					# K(Ψ) oberved
						Kθ_Obs = kunsat.Ψ_2_KUNSAT(option.hydro, Ψ_Obs[iΨ], iZ, hydro)

						Kθ_Log_Obs[iΨ] = log1p(cst.MmS_2_MmH * Kθ_Obs[ipGroup[iZ]])

				end # iΨ
			end # if GroupBool_Select[iZ]
			Of_Kθ = Of_Kθ + (1.0 - stats.NSE_WILMOT(Kθ_Log_Obs , Kθ_Log_Sim)) / sum(GroupBool_Select[iZ])
			end # iZ)
		
	
		return Of_Kθ
		end  # function: OF_KSMODELa
	# --------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SEARCHRANGE
	#		Required by BlackBoxOptim
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SEARCHRANGE(iGroup_Opt, optimKsmodel)
			ParamOpt_Min₂ = copy(optimKsmodel.ParamOpt_Min[iGroup_Opt, 1:optimKsmodel.NparamOpt[iGroup_Opt]])
			ParamOpt_Max₂ = copy(optimKsmodel.ParamOpt_Max[iGroup_Opt, 1:optimKsmodel.NparamOpt[iGroup_Opt]])
		return SearchRange = (collect(zip(Float64.(ParamOpt_Min₂), Float64.(ParamOpt_Max₂))))
		end  # function: SEARCHRANGE
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function X_2_τ(iGroup_Opt, ksmodelτ, optimKsmodel, X)
			for iParam = 1:optimKsmodel.NparamOpt[iGroup_Opt]
				Paramₐ = X[iParam]
				
				# Getting the current values of every layer of the hydro parameter of interest
					vectParam = getfield(ksmodelτ, Symbol(optimKsmodel.ParamOpt[iGroup_Opt, iParam]))

				# Updating the value of the parameters for the layer wanting to optimize by keeping the other values constant
					vectParam[iGroup_Opt] = Paramₐ

				# Putting the updated hydro into ksmodelτ
					setfield!(ksmodelτ, Symbol(optimKsmodel.ParamOpt[iGroup_Opt, iParam]), vectParam)
			end # for loop
		return ksmodelτ
		end  # function: PARAM
	#..................................................................
	
end  # module: optKsModel
# ========================================================================