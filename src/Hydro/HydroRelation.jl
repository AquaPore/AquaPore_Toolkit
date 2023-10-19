# =============================================================
#		MODULE: hydroRealation
# =============================================================
module hydroRelation
import BlackBoxOptim
import ..tool
export σ_2_Ψm, σ_2_θr, FUNCTION_σ_2_Ψm_SOFTWARE

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_Ψm⍰(iZ, hydro)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function σ_2_Ψm(σ₁, Ψσ, Ψm_Min, Ψm_Max; Pσ=3.0)
			Ψm = Ψσ * exp(σ₁ * Pσ)
		return max(min(Ψm , Ψm_Max), Ψm_Min)
		end # function: σ_2_Ψm
	# ----------------------------------------------------------------

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ΨMacMat_FUNC!
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ΨMacMat_FUNC!(hydro₂, iZ; ΨMacMat_Max=100.0, ΨMacMat_Min=10.0, θsMacMat_η_Tresh=0.9)

				θsMacMat_η = min(max((hydro₂.θsMacMat[iZ] - hydro₂.θr[iZ]) / ( hydro₂.θs[iZ] -  hydro₂.θr[iZ]), 0.0), 1.0)

				if θsMacMat_η ≥ θsMacMat_η_Tresh
					return (θsMacMat_η - 1.0) * (ΨMacMat_Max - ΨMacMat_Min) / (θsMacMat_η_Tresh - 1.0) + ΨMacMat_Min
				else
					return ΨMacMat_Max
				end
			end  # function: name
		# ------------------------------------------------------------------

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNCTION_σ_2_Ψm_SOFTWARE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNCTION_σ_2_Ψm_SOFTWARE(hydro₂, iZ, option₂, param; Pσ=3.0)

			# ΨMacMat = ΨMacMat_FUNC!(hydro₂, iZ; ΨMacMat_Max=100.0, ΨMacMat_Min=10.0, θsMacMat_η_Tresh=0.95)

			ΨMacMat = 100.0

			if (option₂.σ_2_Ψm⍰ == "Constrained")
				Ψm_Min = hydroRelation.σ_2_Ψm(hydro₂.σ[iZ], √(ΨMacMat), hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ];Pσ=3)

				Ψm_Max = hydroRelation.σ_2_Ψm(hydro₂.σ[iZ], ΨMacMat, hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ]; Pσ=3)
				
				hydro₂.Ψm[iZ] = tool.norm.∇NORM_2_PARAMETER(hydro₂.Ψm[iZ], Ψm_Min, Ψm_Max)

			elseif (option₂.σ_2_Ψm⍰ == "UniqueRelationship") # <>=<>=<>=<>=<>
				hydro₂.Ψm[iZ] = hydroRelation.σ_2_Ψm(hydro₂.σ[iZ], exp((log(√ΨMacMat) + log(ΨMacMat)) * 0.5), hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ]; Pσ=Pσ)

			end #option.infilt.σ_2_Ψm⍰
		return hydro₂
		end
	# ----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_θr
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function σ_2_θr(hydro₂, iZ; α₁=17.5, α₂=4.0)
		# 	σ_η = (hydro₂.σ[iZ] - hydro₂.σ_Min[iZ]) / (hydro₂.σ_Max[iZ] - hydro₂.σ_Min[iZ]) 	
		# return (hydro₂.θr_Max[iZ] * (1.0 - exp(-α₁ * σ_η ^ α₂))) / (1.0 - exp(-α₁ * 1.0 ^ α₂))
		# end  # function: σ_2_θr

			# α₁_Max =0.75
			# α₁_Min =0.6
			# function σ_2_θr0(hydro₂, iZ; α₁=0.7)
			# 	σ_η = (hydro₂.σ[iZ] - hydro₂.σ_Min[iZ]) / (hydro₂.σ_Max[iZ] - hydro₂.σ_Min[iZ]) 
				
			# 	return hydro₂.θr_Max[iZ] * min(σ_η / α₁, 1.0)
			# end # function: σ_2_θr


			function σ_2_θr(hydro₂, iZ; X1=0.7, X2=2.9)
            Y2 = hydro₂.θr_Max[iZ]
            Y1 = 0.0
            A  = (Y2 - Y1) / (X2 - X1)
            B  = Y1 - X1 * A
			return min(max(A * hydro₂.σ[iZ] + B, 0.0), hydro₂.θr_Max[iZ])
			end # function: σ_2_θr

end  # module: hydroRealation