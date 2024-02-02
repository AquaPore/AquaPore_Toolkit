# =============================================================
#		MODULE: hydroRealation
# =============================================================
module hydroRelation
import BlackBoxOptim
import ..tool
export σ_2_θr, FUNCTION_σ_2_Ψm_SOFTWARE, FUNC_ΨmacMat_2_ΨmMac, FUNC_θsMacMatη_2_ΨmacMat, FUNC_ΨmacMat_2_σMac, FUNC_σ_2_Ψm, FUNC_ΨmMode, FUNC_θsMacMatη


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : name
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function FUNC_θsMacMatη(;θr, θs, θsMacMat)
		return θsMacMat_η = min(max((θsMacMat - θr) / (θs - θr), 0.0), 1.0)
	end  # function: FUNC_θs_θr_2_θsMacMatη
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNC_θsMacMatη_2_ΨmacMat
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function FUNC_θsMacMatη_2_ΨmacMat(;θs, θsMacMat, θr, ΨmacMat_Max=100.0, ΨmacMat_Min=0.0, θsMacMat_η_Tresh=1.0) 

		# 	# θsMacMat_η = FUNC_θsMacMatη(;θr, θs, θsMacMat)

		# 	# if θsMacMat_η ≥ θsMacMat_η_Tresh
		# 	# 	ΨmacMat = ((θsMacMat_η - 1.0) / (θsMacMat_η_Tresh - 1.0)) * (ΨmacMat_Max - ΨmacMat_Min) + ΨmacMat_Min
		# 	# else
		# 	# 	ΨmacMat = ΨmacMat_Max
		# 	# end
			
		# 	ΨmacMat = ΨmacMat_Max

		# 	return ΨmacMat
		# end  # function: FUNC_θsMacMatη_2_ΨmacMat
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNC_ΨmacMat_2_σMac
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNC_ΨmacMat_2_σMac(;ΨmacMat, Pσ_Mac=2)
			# return σMac = log1p(ΨmacMat) / (2.0 * Pσ_Mac)
			return σMac = log(√(ΨmacMat + 1.0)) / Pσ_Mac
		end  # function: FUNC_ΨmacMat_2_σMac
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNC_ΨmacMat_2_ΨmMac
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNC_ΨmacMat_2_ΨmMac(;ΨmacMat, σMac, Option_Mode=false)

			# Option 1: based on Mode of pore size distribution
				# LONG: ΨmMac = exp(log(ΨmacMat) * 0.5 + σMac ^ 2.0)
				# LONG2 ΨmMac = exp(log(sqrt(ΨmacMat)) + σMac ^ 2.0)
				if Option_Mode
					return ΨmMac = √(ΨmacMat + 1.0) * exp(σMac ^ 2.0)
				else 
					return ΨmMac = √(ΨmacMat + 1.0)
				end			
		end  # function: FUNC_ΨmacMat_2_ΨmMac
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNC_σ_2_Ψm
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNC_σ_2_Ψm(;ΨmacMat, σ, Pσ, Ψm_Min=ΨmacMat, Ψm_Max=10.0^8, Option_Mode=false)
			# if Option_Mode
			# 	Ψm = (1.0 + ΨmacMat) * exp(σ * Pσ + σ^2)
			# else
				Ψm = (1.0 + ΨmacMat) * exp(σ * Pσ)
			# end	

			Ψm = max(min(Ψm, Ψm_Max), Ψm_Min) 

		return Ψm
		end # function: FUNC_σ_2_Ψm
	# ----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNCTION_σ_2_Ψm_SOFTWARE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNCTION_σ_2_Ψm_SOFTWARE(hydro₂, iZ, option₂, param; Pσ=3.0, Pσ_Mac=2)


			if option₂.ΨmacMat_2_σMac_ΨmMac
				# ΨmacMat₁ = FUNC_θsMacMatη_2_ΨmacMat(θs=hydro₂.θs[iZ], θsMacMat=hydro₂.θsMacMat[iZ], θr=hydro₂.θr[iZ], ΨmacMat_Max=hydro₂.ΨmacMat[iZ])
				ΨmacMat₁ = hydro₂.ΨmacMat[iZ]
		
				# Deriving σMac
					hydro₂.σMac[iZ] = hydroRelation.FUNC_ΨmacMat_2_σMac(;ΨmacMat=ΨmacMat₁, Pσ_Mac=Pσ_Mac)

				# Deriving ΨmMac
					hydro₂.ΨmMac[iZ] =  hydroRelation.FUNC_ΨmacMat_2_ΨmMac(;ΨmacMat=ΨmacMat₁, σMac=hydro₂.σMac[iZ])
			end

			 ΨmacMat₂ = hydro₂.ΨmacMat[iZ]


			if (option₂.σ_2_Ψm⍰ == "Constrained")
				# Deriving  Ψm 

               Ψm_Min        = hydroRelation.FUNC_σ_2_Ψm(;ΨmacMat=√ΨmacMat₂, σ=hydro₂.σ[iZ], Pσ=Pσ, Ψm_Min=hydro₂.ΨmacMat[iZ], Ψm_Max=hydro₂.Ψm_Max[iZ])

               Ψm_Max        = hydroRelation.FUNC_σ_2_Ψm(;ΨmacMat=ΨmacMat₂, σ=hydro₂.σ[iZ], Pσ=Pσ, Ψm_Min=hydro₂.ΨmacMat[iZ], Ψm_Max=hydro₂.Ψm_Max[iZ])
					
               hydro₂.Ψm[iZ] = tool.norm.∇NORM_2_PARAMETER(hydro₂.Ψm[iZ], Ψm_Min, Ψm_Max)

			elseif (option₂.σ_2_Ψm⍰ == "UniqueRelationship") # <>=<>=<>=<>=<>
				# hydro₂.Ψm[iZ] = hydroRelation.σ_2_Ψm(hydro₂.σ[iZ], exp((log(√ΨmacMat) + log(ΨmacMat)) * 0.5), hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ]; Pσ=Pσ)
				 hydro₂.Ψm[iZ] =  hydroRelation.FUNC_σ_2_Ψm(;ΨmacMat=exp((log(√ΨmacMat₂) + log(ΨmacMat₂)) * 0.5), σ=hydro₂.σ[iZ], Pσ=Pσ, Ψm_Min=hydro₂.Ψm_Min[iZ], Ψm_Max=hydro₂.Ψm_Max[iZ])

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