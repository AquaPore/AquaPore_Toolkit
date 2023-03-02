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
	#		FUNCTION : σ_2_θr
	#		Old
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNCTION_σ_2_Ψm_SOFTWARE(hydro₂, iZ, option₂, param; Pσ=3.0)
			if (option₂.σ_2_Ψm⍰ == "Constrained")
				Ψm_Min = hydroRelation.σ_2_Ψm(hydro₂.σ[iZ], √(param.ΨmacMat), hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ];Pσ=Pσ)

				Ψm_Max = hydroRelation.σ_2_Ψm(hydro₂.σ[iZ], param.ΨmacMat, hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ]; Pσ=Pσ)
				
				hydro₂.Ψm[iZ] = tool.norm.∇NORM_2_PARAMETER(hydro₂.Ψm[iZ], Ψm_Min, Ψm_Max)

			elseif (option₂.σ_2_Ψm⍰ == "UniqueRelationship") # <>=<>=<>=<>=<>
				hydro₂.Ψm[iZ] = hydroRelation.σ_2_Ψm(hydro₂.σ[iZ], (√param.ΨmacMat +  param.ΨmacMat)*0.5, hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ]; Pσ=Pσ)

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


			function σ_2_θr(hydro₂, iZ; X1=0.75, X2=2.9)
            Y2 = hydro₂.θr_Max[iZ]
            Y1 = 0.0
            A  = (Y2 - Y1) / (X2 - X1)
            B  = Y1 - X1 * A
			return max(A * hydro₂.σ[iZ] + B, 0.0)
			end # function: σ_2_θr

			# function σ_2_θr(hydro₂, iZ; X1=0.3, X2=2.6)
         #    Y2 = 1.0
         #    Y1 = 0.0
         #    A  = (Y2 - Y1) / (X2 - X1)
         #    B  = Y1 - X1 * A
			# return hydro₂.θr_Max[iZ] * (max(A * hydro₂.σ[iZ] + B, 0.0)^0.5)
			# end # function: σ_2_θr
	
end  # module: hydroRealation
# ............................................................