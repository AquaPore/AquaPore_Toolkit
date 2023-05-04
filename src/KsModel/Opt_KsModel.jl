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
		function START_OPT_KθMODEL(∑Psd, 🎏_Clay, ClassBool_Select, hydro, ipClass, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param)
				
			# Deriving the feasible range of the τ parameters
				SearchRange = SEARCHRANGE(ipClass, optimKsmodel)

			# Optimisation algorithme, MaxFuncEvals=1000
				Optimization = BlackBoxOptim.bboptimize(X -> OF_KθMODEL(∑Psd, 🎏_Clay, ClassBool_Select, hydro, ipClass, ksmodelτ, NiZ, optim, optimKsmodel, option, param, X; 🎏_IsTopsoil=false, 🎏_RockFragment=false, IsTopsoil=[], RockFragment=[]); SearchRange=SearchRange, NumDimensions=optimKsmodel.NparamOpt[ipClass], TraceMode=:silent)

			# Deriving the optimal τ parameters from X
				X = BlackBoxOptim.best_candidate(Optimization)

			# Putting X parameters into τ
				ksmodelτ = X_2_τ(ipClass, ksmodelτ, optimKsmodel, X)

			# Computing optimal KₛModel
			for iZ=1:NiZ
				if ClassBool_Select[iZ]
					KₛModel[iZ] = θψ_2_KsψModel.KSΨMODEL_START(∑Psd, 🎏_Clay, hydro, ipClass, iZ, ksmodelτ, option, param, 0.0; 🎏_IsTopsoil=false, 🎏_RockFragment=false, IsTopsoil=[], RockFragment=[])
				end
			end
		return KₛModel
		end  # function: START_OPT_KθMODEL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_KθMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_KθMODEL(∑Psd, 🎏_Clay, ClassBool_Select, hydro, ipClass, ksmodelτ, NiZ, optim, optimKsmodel, option, param, X; 🎏_IsTopsoil=false, 🎏_RockFragment=false, IsTopsoil=[], RockFragment=[], KsMinMax=0.005555556)

			# Deriving the optimal τ parameters from X
				ksmodelτ = X_2_τ(ipClass, ksmodelτ, optimKsmodel, X)
				
				Ψ_Obs = param.ksModel.Ψ_Obs		
				N_ΨObs = length(Ψ_Obs)

				Kθ_Log_Obs = fill(0.0::Float64, N_ΨObs)
				Kθ_Log_Sim = fill(0.0::Float64, N_ΨObs)

			# Computing K(Ψ)
			Of_Kθ = 0.0
			for iZ=1:NiZ
			if ClassBool_Select[iZ]
				for iΨ =1:N_ΨObs
					# K(Ψ) simulated
						Kθ_Sim = θψ_2_KsψModel.KSΨMODEL_START(∑Psd, 🎏_Clay, hydro, ipClass, iZ, ksmodelτ, option, param, Ψ_Obs[iΨ]; 🎏_IsTopsoil=false, 🎏_RockFragment=false, IsTopsoil=[], RockFragment=[])

						Kθ_Log_Sim[iΨ] = log(Kθ_Sim)

					# K(Ψ) oberved
						Kθ_Obs = kunsat.Ψ_2_KUNSAT(option.hydro, Ψ_Obs[iΨ], iZ, hydro)
						Kθ_Log_Obs[iΨ] = log(Kθ_Obs)
				end # for iΨ =1:N_ΨObs

				if option.ksModel.Of_KₛModel⍰ == "Wilmot"
					Of_Kθ = Of_Kθ + (1.0 - abs(stats.NSE_WILMOT(Kθ_Log_Obs[1:N_ΨObs], Kθ_Log_Sim[1:N_ΨObs])))
				else
					Of_Kθ = Of_Kθ + stats.RMSE_CONCORDANCE_CORELATION_COEFICIENT(Kθ_Log_Obs[1:N_ΨObs], Kθ_Log_Sim[1:N_ΨObs])
				end
			end # if ClassBool_Select[iZ]
			end # for iZ=1:NiZ		
		return Of_Kθ 
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
	
end  # module: optKsModel
# ========================================================================