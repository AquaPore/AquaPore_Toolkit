# =============================================================
#		module: kunsatModel jesus
# =============================================================
module θψ_2_KsψModel
	import ..cst, ..distribution, ..wrc, ..kunsat
	import QuadGK
	import SpecialFunctions: erfc, erfcinv
	
	export KSΨMODEL_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSΨMODEL_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSΨMODEL_START(∑Psd, 🎏_Clay, hydro, ipClass, iZ, ksmodelτ, option, param, Ψ₁; 🎏_IsTopsoil=false, 🎏_RockFragment=false, IsTopsoil=[], RockFragment=[])

			if 🎏_RockFragment
				RockFragment₁ = RockFragment[iZ]
			else
				RockFragment₁ = 0.0
			end #@isdefined RockFragment

			# if 🎏_IsTopsoil
			# 	IsTopsoil₁ = Int64(IsTopsoil[iZ])
			# else
			# 	IsTopsoil₁ = 1	# Default value				
			# end  # if: @isdefined IsTopsoil

			return KsΨmodel = TORTUOSITYMODELS(∑Psd, 🎏_Clay, hydro, ipClass, iZ, ksmodelτ, option, param, Ψ₁; RockFragment=RockFragment₁, Smap_ImpermClass=[], KsImpClass_Dict=[])
		end  # function: KS_MODEL
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSMODEL_TRADITIONAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSMODEL_TRADITIONAL(Se, T1, T1Mac, T2, T2Mac, T3, T3Mac, θr, θs, θsMacMat, σ, σMac, Ψm, ΨmMac)

			Kunsat_Mat = T1 * ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

			Kunsat_Mac = T1Mac * ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 
		return KsModel = Kunsat_Mat + Kunsat_Mac
		end  # function: KS_MODEL
	# ------------------------------------------------------------------
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KsΨMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KsΨMODEL(hydro, iZ::Int64, optionₘ, T1, T1Mac, T2, T2Mac, T3, T3Mac, θr, θs, θsMacMat, σ, σMac, Ψ₁::Float64, Ψm, ΨmMac)

			# Matrix ====	
				θ_Mat = 0.5 * (θsMacMat - θr) * erfc((log( Ψ₁ / Ψm)) / (σ * √2.0)) + θr
				
				Se_Mat = (θ_Mat- θr) / (θs - θr)
				
				Ks_Mat = T1 * cst.KunsatModel * π * ((θsMacMat - θr) * ((cst.Y / Ψm) ^ T2) * exp(((T2 * σ) ^ 2.0) / 2.0)) ^ T3
				
				Kunsat_Mat = Ks_Mat * √Se_Mat * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0
			
			# Macropore ===
				θ_Mac = 0.5 * (θs - θsMacMat) * erfc((log(Ψ₁ / ΨmMac)) / (σMac * √2.0))
				
				Se_Mac = θ_Mac / (θs - θr)

				Ks_Mac = T1Mac * cst.KunsatModel * π * ((θs - θsMacMat) * ((cst.Y / ΨmMac) ^ T2Mac) * exp(((T2Mac * σMac) ^ 2.0) / 2.0)) ^ T3Mac

				Kunsat_Mac = Ks_Mac * √Se_Mac * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0
		return K_Ψ = Kunsat_Mat + Kunsat_Mac
		end  # function: KsΨMODEL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KsΨMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KsΨMODEL_CLAY(hydro, iZ::Int64, optionₘ, T1, T1Mac, T2, T2Mac, T3, T3Mac, Tclay, θr, θs, θsMacMat, σ, σMac, Ψ₁::Float64, Ψm, ΨmMac)

			# Se ===
				Se = wrc.kg.Ψ_2_SeDual(optionₘ, Ψ₁, iZ, hydro)

			# Matrix ====	
				Ks_Mat = T1 * cst.KunsatModel * π * (((θsMacMat - θr) ^ (Tclay / T3)) * ((cst.Y / Ψm) ^ T2) * exp(((T2 * σ) ^ 2.0) / 2.0)) ^ T3

				Kunsat_Mat = Ks_Mat * √Se * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0

			# Macropore ===
				Ks_Mac = T1Mac * cst.KunsatModel * π * ((θs - θsMacMat)  * ((cst.Y / ΨmMac) ^ T2Mac) * exp(((T2Mac * σMac) ^ 2.0) / 2.0)) ^ T3Mac

				Kunsat_Mac = Ks_Mac * √Se * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0

	return K_Ψ = Kunsat_Mat + Kunsat_Mac
	end  # function: KsΨMODEL
	# ------------------------------------------------------------------

	# ^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^__^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^__^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TORTUOSITYMODELS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TORTUOSITYMODELS(∑Psd, 🎏_Clay::Bool, hydro, ipClass, iZ::Int64, ksmodelτ, option, param, Ψ₁; RockFragment=0.0, θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψm=hydro.Ψm[iZ], σ=hydro.σ[iZ], θsMacMat=hydro.θsMacMat[iZ], ΨmMac=hydro.ΨmMac[iZ], σMac=hydro.σMac[iZ], τ₁ₐ=ksmodelτ.τ₁ₐ[ipClass],τclay₀=ksmodelτ.τclay₀[ipClass], τ₂ₐ=ksmodelτ.τ₂ₐ[ipClass], τclayₘₐₓ=ksmodelτ.τclayₘₐₓ[ipClass], τ₃ₐ=ksmodelτ.τ₃ₐ[ipClass], τclayΔθsr=ksmodelτ.τclayΔθsr[ipClass], τ₁ₐMac=ksmodelτ.τ₁ₐMac[ipClass],τclay₀Mac=ksmodelτ.τclay₀Mac[ipClass], τ₂ₐMac=ksmodelτ.τ₂ₐMac[ipClass], τclayₘₐₓMac=ksmodelτ.τclayₘₐₓMac[ipClass], τ₃ₐMac=ksmodelτ.τ₃ₐMac[ipClass], τclayΔθsrMac=ksmodelτ.τclayΔθsrMac, RockFragment_Treshold=0.4, Smap_ImpermClass=[], KsImpClass_Dict=[] )

			# Determine when Ks increases for increasing RockFragment	
				if RockFragment > RockFragment_Treshold
					θr, θs, θsMacMat = ROCKCORRECTION(RockFragment, RockFragment_Treshold, θr, θs, θsMacMat)
				end

			# MODEL 0 ====
			# Original model	
			if option.ksModel.KₛModel⍰=="KsΨmodel_0" # ===
				T2_Max = 2.0

				# Transforming matrix
					T1 =  10.0 ^ (τ₁ₐ / (τ₁ₐ - 1.0))
					T2 = T2_Max * (1.0 - τ₂ₐ)
					T3 = 1.0 / (1.0 - τ₃ₐ)

				# Transforming macro
					T1Mac = 10.0 ^ (τ₁ₐMac / (1.0 - τ₁ₐMac))
					T2Mac = T2_Max * (1.0 - τ₂ₐMac)
					T3Mac = 1.0 / (1.0 - τ₃ₐMac)	
		
				# Ks model
					KₛModel = cst.KunsatModel * QuadGK.quadgk(Se -> KSMODEL_TRADITIONAL(Se, T1, T1Mac, T2, T2Mac, T3, T3Mac, θr, θs, θsMacMat, σ, σMac, Ψm, ΨmMac), 0.0, 0.9999; rtol=1.0E-3)[1]

					return Kunsat = KₛModel * kunsat.kg.Ψ_2_KUNSAT(option.hydro, Ψ₁, iZ, hydro) / hydro.Ks[iZ]
				
			# MODEL 2 ====	
			# Rekationship between macro and matrix		
			elseif option.ksModel.KₛModel⍰=="KsΨmodel_1" # ===	
				T2_Max = 3.0; T3_Max = 4.0
				# MATRIX 
					# Tortuosity T1
						T1 = 10.0 ^ (τ₁ₐ / (τ₁ₐ - 1.0))
						
					# Tortuosity T2
						T2 = T2_Max * (1.0 - τ₂ₐ)

					# Tortuosity T3	
						T3 = T3_Max * (1.0 - τ₃ₐ)			

				# MACRO
					# Tortuosity T1Mac
						T1Mac = 10.0 ^ (τ₁ₐMac / (τ₁ₐMac - 1.0))
					
					# Tortuosity T2Mac
						T2Mac = T2_Max * (1.0 - τ₂ₐMac)
					
					# Tortuosity T3Mac
						T3Mac = T3_Max * (1.0 - τ₃ₐMac)			 											
			return KsΨMODEL(hydro, iZ, option.hydro, T1, T1Mac, T2, T2Mac, T3, T3Mac, θr, θs, θsMacMat, σ, σMac, Ψ₁, Ψm, ΨmMac)
				
	
			# MODEL 2 ====	
			# Rekationship between macro and matrix
			elseif option.ksModel.KₛModel⍰=="KsΨmodel_2"
				# CLAY FUNCTION
					# CLAY FUNCTION					
					Tclay = TORTUOSITY_CLAY(∑Psd, 🎏_Clay, hydro, iZ, option, param, τclay₀, τclayₘₐₓ, τclayΔθsr)
				
				# MATRIX 
						T2_Max = 3.0; T3_Max = 4.0
					# Tortuosity T1
						T1 =  10.0 ^ (τ₁ₐ / (τ₁ₐ - 1.0))
						# T1 = 10.0 ^ - (1.0 / (Tclay₀ * τ₁ₐ) -1.0)
						
					# Tortuosity T2
						T2 = T2_Max * (1.0 - τ₂ₐ)

					# Tortuosity T3	
						T3 = T3_Max * (1.0 - τ₃ₐ)

				# MACRO
					# Tortuosity T1Mac
						T1Mac = 10.0 ^ (τ₁ₐMac / (τ₁ₐMac - 1.0))
					
					# Tortuosity T2Mac
						T2Mac = T2_Max * (1.0 - τ₂ₐMac)
					
					# Tortuosity T3Mac
						T3Mac = T3_Max * (1.0 - τ₃ₐMac)
							
				return KsΨMODEL_CLAY(hydro, iZ, option.hydro, T1, T1Mac, T2, T2Mac, T3, T3Mac, Tclay, θr, θs, θsMacMat, σ, σMac, Ψ₁, Ψm, ΨmMac)	
		
			
			# MODEL 3 ====	
			# Rekationship between macro and matrix
			elseif option.ksModel.KₛModel⍰=="KsΨmodel_3"
				# CLAY FUNCTION					
					Tclay = TORTUOSITY_CLAY(∑Psd, 🎏_Clay, hydro, iZ, option, param, τclay₀, τclayₘₐₓ, τclayΔθsr)
						
				# MATRIX ----------------------------
						T2_Max = 3.0; T3_Max = 4.0
					# Tortuosity T1
						T1 =  10.0 ^ (τ₁ₐ / (τ₁ₐ - 1.0))
						
					# Tortuosity T2
						T2 = T2_Max * (1.0 - τ₂ₐ)

					# Tortuosity T3	
						T3 = T3_Max * (1.0 - τ₃ₐ)

				# MACRO ----------------------------
					# Tortuosity T1Mac
						T1Mac = 10.0 ^ (τ₁ₐ / (τ₁ₐ - 1.0))
					
					# Tortuosity T2Mac
						T2Mac = T2_Max * (1.0 - τ₂ₐMac)
					
					# Tortuosity T3Mac
						T3Mac = T3_Max * (1.0 - τ₃ₐMac)
						
				return KsΨMODEL_CLAY(hydro, iZ, option.hydro, T1, T1Mac, T2, T2Mac, T3, T3Mac, Tclay, θr, θs, θsMacMat, σ, σMac, Ψ₁, Ψm, ΨmMac)	
			
			# Model 3Unimodal	
			elseif option.ksModel.KₛModel⍰=="KsΨmodel_3Unimodal"
				# CLAY FUNCTION
					# CLAY FUNCTION					
					Tclay = TORTUOSITY_CLAY(∑Psd, 🎏_Clay, hydro, iZ, option, param, τclay₀, τclayₘₐₓ, τclayΔθsr)
						
				# MATRIX 
						T2_Max = 3.0; T3_Max = 3.0
					# Tortuosity T1
						T1 = 10.0 ^ (τ₁ₐ / (τ₁ₐ - 1.0))

						# T1 = 10.0 ^ - (1.0 / (Tclay₀ * τ₁ₐ) -1.0)
						
					# Tortuosity T2
						T2 = T2_Max * (1.0 - τ₂ₐ)

					# Tortuosity T3	
						T3 = T3_Max * (1.0 - τ₃ₐ)

				# MACRO
					# Tortuosity T1Mac
						τ₁ₐMac = 0.0
						T1Mac = 10.0 ^ (τ₁ₐMac / (τ₁ₐMac - 1.0))
					
					# Tortuosity T2Mac
						τ₂ₐMac = 0.0
						T2Mac = T2_Max * (1.0 - τ₂ₐMac) 
					
					# Tortuosity T3Mac
						τ₃ₐMac = 0.0
						T3Mac = T3_Max * (1.0 - τ₃ₐMac)
							
				return KsΨMODEL_CLAY(hydro, iZ, option.hydro, T1, T1Mac, T2, T2Mac, T3, T3Mac, Tclay, θr, θs, θsMacMat, σ, σMac, Ψ₁, Ψm, ΨmMac)	
			
			elseif option.ksModel.KₛModel⍰=="KsΨmodel_4" # ===
				# Transformation matrix

				# CLAY MODEL 
					# Reducing with σ
						σclay = 2.3 # 2.3

						X_σ₁ = 0
						Y_σ₁ = 1.0 
						X_σ₂ = 3.7 - σclay
						Y_σ₂ = τclay₀
						α  = (Y_σ₂ - Y_σ₁) / (X_σ₂ - X_σ₁)
						B  = Y_σ₁ - X_σ₁ * α 

					Tσ = max(min(α * (σ - σclay) + B, 1.0), Y_σ₂)

					# Reducing with θs - θr
						X_θs₁ = 0.0
						Y_θs₁ = 1.0
						X_θs₂ = 0.6
						Y_θs₂ = 0.0
						α = (Y_θs₂ - Y_θs₁) / (X_θs₂  - X_θs₁)
						Β  = Y_θs₁ - X_θs₁ * α

						TθsMacMat = max(min(α * (θsMacMat - θr) + Β, 1.0), Y_θs₂)
					
					# T1 TORTUSOSITY MODEL
						T1 = 10.0 ^ (τ₁ₐ / (τ₁ₐ - 1.0))
						if σ ≥ σclay 
							if (θsMacMat - θr) ≥ 0.25
								T1 = T1 * TθsMacMat ^ τ₂ₐ
							end
						end
			
				# Tortuosity T2
					T2_Min = 1.0; T2_Max = 3.0
					T2 = ((T2_Min - T2_Max) * τ₂ₐ + T2_Max)

					if σ ≥ σclay 					
						T2 = T2 * Tσ 
					end

				# Tortuosity T3
					T3 = τ₃ₐ

			# Transformation macro
				T1Mac = 10.0 ^ (τ₁ₐMac / (τ₁ₐMac - 1.0))
				T2Mac = (T2_Min - T2_Max) * τ₂ₐMac + T2_Max
				T3Mac = τ₃ₐMac							
			
			return KsΨMODEL(hydro, iZ, option.hydro, T1, T1Mac, T2, T2Mac, T3, T3Mac, θr, θs, θsMacMat, σ, σMac, Ψ₁, Ψm, ΨmMac)
				
			else
				error("option.ksModel.KₛModel⍰ = $(option.ksModel.KₛModel⍰) is not yet implemented try <KsModel_Traditional>; <KsModel_Tσ>; <KsModel_New>; <KsModel_NewSimplified> ")
				
			end  # if: Model=="Model?"
		end  # function: TORTUOSITYMODELS 


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ROCKCORRECTION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ROCKCORRECTION(RockFragment, RockFragment_Treshold, θr, θs, θsMacMat)
				RockFragment2 = max(2.0 * RockFragment_Treshold - RockFragment, 0.0)

				θs = (θs / (1.0 - RockFragment)) * (1.0 - RockFragment2)
				
				θsMacMat = (θsMacMat / (1.0 - RockFragment)) * (1.0 - RockFragment2)

				θr = (θr / (1.0 - RockFragment)) * (1.0 - RockFragment2)		
		return θr, θs, θsMacMat
		end  # function: ROCKCORRECTION
	# ------------------------------------------------------------------

	# =====================================================================================================================
	# =====================================================================================================================

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TORTUOSITY_CLAY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TORTUOSITY_CLAY(∑Psd, 🎏_Clay, hydro, iZ, option, param, τclay₀, τclayₘₐₓ, τclayΔθsr)

			# If we have clay information derived from PSD
			if 🎏_Clay
				Clay = ∑Psd[iZ,1]
			
			# Rough modelling % of clay [0-1]
			# Correlation between clay particle and Ψ
			else
				# Ψ_Clay =  160000.0 * ( ( (cst.Y  / 0.002) - (cst.Y / 0.5) ) / ((cst.Y  / 0.001) - (cst.Y  / 0.5)) ) ^ 2.0
				Ψ_Clay = param.psd.imp.Ψ_Max * (((cst.Y / 0.002) - (cst.Y / 0.5) ) / ((cst.Y / 0.002) - (cst.Y / 0.5))) ^ param.psd.imp.λ 

				Clay = wrc.Ψ_2_SeDual(option.hydro, Ψ_Clay, iZ, hydro)
			end # 🎏_Clay
			
			X_Clay₁ =  τclay₀

			Clayₙ = max(Clay - X_Clay₁, 0.0) / (1.0 - X_Clay₁)

			ΔθsMacθr = hydro.θsMacMat[iZ] - hydro.θr[iZ]

			ΔθsMacθrₙ =  max(ΔθsMacθr - τclayΔθsr , 0.0) / (1.0 - τclayΔθsr)

			Tclay_Max =  1.0 + ΔθsMacθrₙ * (τclayₘₐₓ - 1.0) 

		return Tclay = Tclay_Max - (Tclay_Max - 1.0) * cos(Clayₙ * π * 0.5) 
		end				
	# ------------------------------------------------------------------
end  # module θψ_2_KsψModel
