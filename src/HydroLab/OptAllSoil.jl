# =============================================================
#		module: optAllSoil
# =============================================================
module optAllSoil
import ..ofHydrolab, ..optimize, ..optIndivSoil, ..ofHydrolab
using BlackBoxOptim
export OPTIMIZE_ALLSOILS

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPTIMIZE_ALLSOILS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function OPTIMIZE_ALLSOILS(;∑Psd::Vector{Any}, hydro::Main.hydroStruct.KOSUGI, hydroOther::Main.hydroStruct.HYDRO_OTHER, K_KΨobs::Matrix{Float64}, N_KΨobs=1, N_θΨobs::Vector{Int64}, NiZ::Int64, optim::Main.reading.OPTIM, optimAllSoils::Main.reading.OPTIM, option::Main.options.OPTION, optionₘ::Main.options.HYDRO, param::Main.params.PARAM, θ_θΨobs::Matrix{Float64}, θϵ=0.005::Float64, Ψ_KΨobs::Matrix{Float64}, Ψ_θΨobs::Matrix{Float64})


# 		const MyFitnessGoal = 30.0

# function myfitnessgoalachieved(oc)
#     best_fitness(oc) < MyFitnessGoal
# end

# function cbearlystopping(oc)
#     if myfitnessgoalachieved(oc)
#         BlackBoxOptim.shutdown!(oc)
#     end
# end

		SearchRange_AllSoils = optimize.SEARCHRANGE(optionₘ, optimAllSoils)

# optimAllSoils.InitialGuess

		Optimization = BlackBoxOptim.bboptimize(X -> optAllSoil.OF_HYDROLAB(;∑Psd, hydro, hydroOther, K_KΨobs, N_KΨobs, N_θΨobs, NiZ, optim, optimAllSoils, option, optionₘ, param, X, θ_θΨobs, θϵ, Ψ_KΨobs, Ψ_θΨobs) ; SearchRange=SearchRange_AllSoils, NumDimensions=optimAllSoils.NparamOpt, TraceMode=:silent)
		# MaxTime=10, 

		# Best parameter set
			X = BlackBoxOptim.best_candidate(Optimization)
			hydro = optAllSoil.PARAM_2_hydro(hydro, NiZ, optimAllSoils, optionₘ, param, X)
		 
	return hydro
	end  # function: OPTIMIZE_ALLSOILS
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_HYDROLAB
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_HYDROLAB(;∑Psd::Vector{Any}, hydro::Main.hydroStruct.KOSUGI, hydroOther::Main.hydroStruct.HYDRO_OTHER, K_KΨobs::Matrix{Float64}, N_KΨobs=1, N_θΨobs::Vector{Int64}, NiZ::Int64, optim::Main.reading.OPTIM, optimAllSoils::Main.reading.OPTIM, option::Main.options.OPTION, optionₘ::Main.options.HYDRO, param::Main.params.PARAM, X, θ_θΨobs::Matrix{Float64}, θϵ=0.005::Float64, Ψ_KΨobs::Matrix{Float64}, Ψ_θΨobs::Matrix{Float64})

		# New optimized which are put into the matching veg or hydro parameters
			hydro = optAllSoil.PARAM_2_hydro(hydro, NiZ, optimAllSoils, optionₘ, param, X)

			hydro, hydroOther, Of_Sample = optIndivSoil.OPTIMIZE_INDIVIDUALSOILS(;∑Psd, hydro, hydroOther, K_KΨobs, N_KΨobs, N_θΨobs, NiZ, optim, option, optionₘ, param, θ_θΨobs, θϵ=0.005, Ψ_KΨobs, Ψ_θΨobs)

			OF = ofHydrolab.OF_ALLSOILS(NiZ, Of_Sample)

		println(" Of =  ", trunc(OF,digits=3), "\n")
		return OF
		end  # function: name
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PARAM_2_hydro
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PARAM_2_hydro(hydro, NiZ, optimAllSoils, optionₘ, param, X)

			for iParam = 1:optimAllSoils.NparamOpt	
				# Determening if parameters are Log transformed
					if (optimAllSoils.ParamOpt_LogTransform[iParam]) && !(optimAllSoils.ParamOpt[iParam]=="Ψm" && optionₘ.σ_2_Ψm⍰ == "Constrained")
						Paramₐ = expm1(X[iParam])
					else
						Paramₐ = X[iParam]
					end  # if: optim.ParamOpt_LogTransform

				# Getting the current values of every layer of the hydro parameter of interest
					vectParam = getfield(hydro, Symbol(optimAllSoils.ParamOpt[iParam]))

				# Updating the value of the parameters for the soil wanting to optimize by keeping the values constant
					for iZ= 1:NiZ
						vectParam[iZ] = Paramₐ
					end

				# Putting the updated hydro into hydro
					setfield!(hydro, Symbol(optimAllSoils.ParamOpt[iParam]), vectParam)
			end # for loop

		return hydro
		end  # function: PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
end  # module optAllSoil
# ............................................................