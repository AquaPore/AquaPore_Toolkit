# =============================================================
#		MODULE: hydroRealation
# =============================================================
module hydroRelation
import BlackBoxOptim
import ..tool
export σ_2_θr, FUNCTION_σ_2_Ψm_SOFTWARE, FUNC_ΨMacMat_2_ΨmMac, FUNC_θsMacMatη_2_ΨMacMat, FUNC_ΨMacMat_2_σMac, FUNC_σ_2_Ψm, FUNC_ΨmMode


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNC_θsMacMatη_2_ΨMacMat
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNC_θsMacMatη_2_ΨMacMat(;θs, θsMacMat, θr, ΨMacMat_Max=70.0, ΨMacMat_Min=0.0, θsMacMat_η_Tresh=0.95) 

			θsMacMat_η = min(max((θsMacMat - θr) / (θs - θr), 0.0), 1.0)

			if θsMacMat_η > θsMacMat_η_Tresh
				ΨMacMat = ((θsMacMat_η - 1.0) / (θsMacMat_η_Tresh - 1.0)) * (ΨMacMat_Max - ΨMacMat_Min) + ΨMacMat_Min
			else
				ΨMacMat = ΨMacMat_Max
			end
		return ΨMacMat
		end  # function: FUNC_θsMacMatη_2_ΨMacMat
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNC_ΨMacMat_2_σMac
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNC_ΨMacMat_2_σMac(;ΨMacMat, Pσ_Mac=2)
			# σMac = log(√ΨMacMat) / Pσ
			return log1p(ΨMacMat) / (2.0 * Pσ_Mac)
		end  # function: FUNC_ΨMacMat_2_σMac
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNC_ΨMacMat_2_ΨmMac
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNC_ΨMacMat_2_ΨmMac(;ΨMacMat, σMac)

			# Option 1: based on Mode of pore size distribution
				# LONG: ΨmMac = exp(log(ΨMacMat) * 0.5 + σMac ^ 2.0)
				# LONG2 ΨmMac = exp(log(sqrt(ΨMacMat)) + σMac ^ 2.0)

				# return √ΨMacMat * exp(σMac ^ 2.0)

			# Option 2:
				return √ΨMacMat						
		end  # function: FUNC_ΨMacMat_2_ΨmMac
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNC_σ_2_Ψm
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNC_σ_2_Ψm(;ΨMacMat, σ, Pσ, Ψm_Min=ΨMacMat, Ψm_Max=10.0^8)

			# Option 1 based on Mode of pore size distribution
				# Ψm = (1.0 + ΨMacMat) * exp(σ * Pσ + σ^2)

			# OPtion 2 based on θ(ψ)
				Ψm = (1.0 + ΨMacMat) * exp(σ * Pσ)

				# Ψm = min(Ψm, Ψm_Max)
		return Ψm
		end # function: FUNC_σ_2_Ψm
	# ----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNCTION_σ_2_Ψm_SOFTWARE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNCTION_σ_2_Ψm_SOFTWARE(hydro₂, iZ, option₂, param; Pσ=3.0)

			ΨMacMat = FUNC_θsMacMatη_2_ΨMacMat(θs=hydro₂.θs[iZ], θsMacMat=hydro₂.θsMacMat[iZ], θr=hydro₂.θr[iZ])
	
			# Deriving σMac
				hydro₂.σMac[iZ] = hydroRelation.FUNC_ΨMacMat_2_σMac(;ΨMacMat, Pσ_Mac=2)

			# Deriving ΨmMac
				hydro₂.ΨmMac[iZ] =  hydroRelation.FUNC_ΨMacMat_2_ΨmMac(;ΨMacMat, σMac=hydro₂.σMac[iZ])

			if (option₂.σ_2_Ψm⍰ == "Constrained")
				# Deriving  Ψm 
               ΨMacMat₂      = 70.0

               Ψm_Min        = hydroRelation.FUNC_σ_2_Ψm(;ΨMacMat=√ΨMacMat₂, σ=hydro₂.σ[iZ], Pσ=3.0, Ψm_Min=hydro₂.Ψm_Min[iZ], Ψm_Max=hydro₂.Ψm_Max[iZ])

               Ψm_Max        = hydroRelation.FUNC_σ_2_Ψm(;ΨMacMat=ΨMacMat₂, σ=hydro₂.σ[iZ], Pσ=3.0, Ψm_Min=hydro₂.Ψm_Min[iZ], Ψm_Max=hydro₂.Ψm_Max[iZ])
					
               hydro₂.Ψm[iZ] = tool.norm.∇NORM_2_PARAMETER(hydro₂.Ψm[iZ], Ψm_Min, Ψm_Max)


			elseif (option₂.σ_2_Ψm⍰ == "UniqueRelationship") # <>=<>=<>=<>=<>
				# hydro₂.Ψm[iZ] = hydroRelation.σ_2_Ψm(hydro₂.σ[iZ], exp((log(√ΨMacMat) + log(ΨMacMat)) * 0.5), hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ]; Pσ=Pσ)
				 hydro₂.Ψm[iZ] =  hydroRelation.FUNC_σ_2_Ψm(;ΨMacMat=exp((log(√ΨMacMat) + log(ΨMacMat)) * 0.5), σ=hydro₂.σ[iZ], Pσ=Pσ, Ψm_Min=hydro₂.Ψm_Min[iZ], Ψm_Max=hydro₂.Ψm_Max[iZ])

			end #option.infilt.σ_2_Ψm⍰
		return hydro₂
		end
	# ----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNC_ΨmMode
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNC_ΨmMode(;Ψm₀ , σ₀)
		# exp(log(Ψm₀) - σ₀^2)
			Ψm_Mode = Ψm₀ * exp(- σ₀^2.0)	
		return Ψm_Mode
		end  # function: FUNC_ΨmMode
	# ------------------------------------------------------------------


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

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_θr
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function σ_2_θr(hydro₂, iZ; X1=0.7, X2=2.9)
			Y2 = hydro₂.θr_Max[iZ]
			Y1 = 0.0
			A  = (Y2 - Y1) / (X2 - X1)
			B  = Y1 - X1 * A
		return min(max(A * hydro₂.σ[iZ] + B, 0.0), hydro₂.θr_Max[iZ])
		end # function: σ_2_θr
	# ----------------------------------------------------------------	

end  # module: hydroRealation