module ofHydrolab
	import..stats, ..wrc, ..kunsat
	export  OF_WRC_KUNSAT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PENALTY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function PENALTY_KUNSAT(hydro, iZ, optionₘ; Macro_Perc_Max=0.4)

		# 	if optionₘ.HydroModel⍰ =="Kosugi" && optionₘ.KosugiModel_KΨ⍰=="ΨmacMat"
		# 		# COMPUTING MACROPORE %
		# 			Ta, Tb, Tc, TaMac, TbMac, TcMac = kunsat.kg.TORTUOSITY(; σ=hydro.σ[iZ], τa=hydro.τa[iZ], τaMac=hydro.τaMac[iZ], τb=hydro.τb[iZ], τbMac=hydro.τbMac[iZ], τc=hydro.τc[iZ], τcMac=hydro.τcMac[iZ])

		# 			KsMac, KsMat= kunsat.kg.KS_MATMAC_ΨmacMat(hydro.θs[iZ], hydro.θsMacMat[iZ], hydro.θr[iZ], hydro.Ψm[iZ], hydro.σ[iZ], σMac=hydro.σMac[iZ],hydro.ΨmMac[iZ], hydro.σMac[iZ], hydro.Ks[iZ], Tb, Tc, TbMac, TcMac, optionₘ.KosugiModel_KΨ⍰)

		# 			Macro_Perc = KsMac / hydro.Ks[iZ]

		# 			Penalty_Ks = max(Macro_Perc - Macro_Perc_Max, 0.0)

		# 			return Penalty_Ks
		# 	else
		# 			return Penalty_Ks = 0.0
		# 	end
		# end  # function: PENALTY
	#--------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : WEIGHT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function WEIGHTING(; iZ, N_KΨobs, N_θΨobs, Wof_Max=0.8, Wof_Min=0.2, Ψ_KΨobs, Ψ_θΨobs)

			Ψ_θΨobs_Min = log10(minimum(abs.(Ψ_θΨobs[iZ, 1:N_θΨobs[iZ]])) + 1.0)
			Ψ_θΨobs_Max = log10(maximum(abs.(Ψ_θΨobs[iZ, 1:N_θΨobs[iZ]])) + 1.0)	

			Ψ_KΨobs_Min = log10(minimum(abs.(Ψ_KΨobs[iZ, 1:N_KΨobs[iZ]])) + 1.0)
			Ψ_KΨobs_Max = log10(maximum(abs.(Ψ_KΨobs[iZ, 1:N_KΨobs[iZ]])) + 1.0)
			
			Wof = 0.4 * (Ψ_θΨobs_Max - Ψ_θΨobs_Min + 1.0) / (Ψ_KΨobs_Max - Ψ_KΨobs_Min + 1.0)
			Wof = max(min(Wof, Wof_Max), Wof_Min)

		return Wof
		end  # function: WEIGHTING
	# ------------------------------------------------------------------

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_WRC_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
		function OF_WRC_KUNSAT(hydro, iZ::Int64, N_θΨobs::Vector{Int64}, Of_Sample::Vector{Float64}, optim, optionₘ, θ_θΨobs::Matrix{Float64}, Ψ_θΨobs::Matrix{Float64}; K_KΨobs=[0.0 0.0; 0.0 0.0]::Matrix{Float64}, N_KΨobs=[0]::Vector{Int64}, Ψ_KΨobs=[0.0 0.0; 0.0 0.0]::Matrix{Float64}, Wof_Min=0.1::Float64, Wof_Max=0.9::Float64)
			
			# === OF θΨ ====
				θ_Obs = fill(0.0::Float64, N_θΨobs[iZ])
				θ_Sim = fill(0.0::Float64, N_θΨobs[iZ])

				# IF we do not optimise θs than we do not use the first obeserved point which is at Ψ=0
					if "θs" ∉ optim.ParamOpt
						iStart = 2
					else
						iStart = 1
					end

				for iΨ = 1:N_θΨobs[iZ]
					θ_Obs[iΨ] = θ_θΨobs[iZ,iΨ]
					θ_Sim[iΨ] = wrc.Ψ_2_θ(optionₘ, Ψ_θΨobs[iZ,iΨ], iZ, hydro)
				end # for iΨ = 1:N_θΨobs[iZ]

				Of_θΨ = stats.NSE_MINIMIZE(θ_Obs[iStart:N_θΨobs[iZ]], θ_Sim[iStart:N_θΨobs[iZ]])

			# === OF Kunsat ====
			if "Ks" ∈ optim.ParamOpt
				# Weighting algorithm
					Wof =  WEIGHTING(; iZ, N_KΨobs, N_θΨobs, Wof_Max=Wof_Max, Wof_Min=Wof_Min, Ψ_KΨobs, Ψ_θΨobs)
					
				Kunsat_Obs_Ln = fill(0.0::Float64, N_KΨobs[iZ])
				Kunsat_Sim_Ln = fill(0.0::Float64, N_KΨobs[iZ])

				for iΨ = 1:N_KΨobs[iZ]
					Kunsat_Obs_Ln[iΨ] = log10(K_KΨobs[iZ,iΨ])
						
					Kunsat_Sim_Ln[iΨ] = log10(kunsat.KUNSAT_θΨSe(optionₘ, Ψ_KΨobs[iZ,iΨ], iZ, hydro))
				end # for iΨ = 1:N_KΨobs[iZ]

				Of_Kunsat = stats.NSE_MINIMIZE(Kunsat_Obs_Ln[1:N_KΨobs[iZ]], Kunsat_Sim_Ln[1:N_KΨobs[iZ]])

				Of_Sample[iZ] = Wof * Of_θΨ + (1.0 - Wof) * Of_Kunsat

			else		
				Of_Sample[iZ] = Of_θΨ
				Of_Kunsat = 0.0
			end #  "Ks" ∈ optim.ParamOpt
		return Of_Sample, Of_θΨ, Of_Kunsat
		end # function OF_WRC_KUNSAT
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_ALLSOILS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_ALLSOILS(NiZ::Int64, Of_Sample::Vector{Float64})
			Of_AllSoil = 0.0::Float64
			for iZ=1:NiZ		
				Of_AllSoil += round(max(min(Of_Sample[iZ], 1.0), 0.0), digits=4)
			end
		return Of_AllSoil = Of_AllSoil / NiZ
		end  # function: OF_ALL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_RMSE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
		function OF_RMSE(option, optionₘ, iZ, θ_θΨobs, Ψ_θΨobs, N_θΨobs, K_KΨobs, Ψ_KΨobs, N_KΨobs, hydro, optim) 

		# === OF θΨ ====
			θ_Sim = fill(0.0::Float64, N_θΨobs[iZ])

			for iΨ = 1:N_θΨobs[iZ]
				θ_Sim[iΨ] = wrc.Ψ_2_θ(optionₘ, Ψ_θΨobs[iZ,iΨ], iZ, hydro)
			end # for iΨ = 1:N_θΨobs[iZ]

			Rmse_θΨ = stats.RMSE(θ_θΨobs[iZ, 1:N_θΨobs[iZ]], θ_Sim[1:N_θΨobs[iZ]])

		# === OF Kunsat ====
			if "Ks" ∈ optim.ParamOpt ||option.run.HydroLabθΨ⍰=="Run"
				if  "Ks" ∈ optim.ParamOpt
					iStart = 1
				else
					iStart = 2
				end

				Kunsat_Obs_Ln = fill(0.0::Float64, N_KΨobs[iZ])
				Kunsat_Sim_Ln = fill(0.0::Float64, N_KΨobs[iZ])

				for iΨ = iStart:N_KΨobs[iZ]
					Kunsat_Obs_Ln[iΨ] = log10(K_KΨobs[iZ,iΨ])
						
					Kunsat_Sim_Ln[iΨ] = log10(kunsat.KUNSAT_θΨSe(optionₘ, Ψ_KΨobs[iZ,iΨ], iZ, hydro))
				end # for iΨ = 1:N_KΨobs[iZ]

				Rmse_KΨ = stats.RMSE(Kunsat_Obs_Ln[iStart:N_KΨobs[iZ]], Kunsat_Sim_Ln[iStart:N_KΨobs[iZ]])

				Rmse = (Rmse_θΨ + Rmse_KΨ) * 0.5
			else		
				Rmse = Rmse_θΨ
				Rmse_KΨ = 0.0
			end #  "Ks" ∈ optim.ParamOpt

	return Rmse, Rmse_KΨ, Rmse_θΨ
	end # OF_RMSE

end # module ofHydaulic
# ............................................................