# =============================================================
#		module: kunsatModel jesus
# =============================================================
module θψ_2_KsψModel
	import ..cst, ..distribution, ..wrc, ..kunsat
	import QuadGK, Polynomials
	import SpecialFunctions: erfc, erfcinv
	
	export KSΨMODEL_START, KsΨMODEL_NOINTEGRAL, TORTUOSITY_CLAY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSΨMODEL_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSΨMODEL_START(∑Psd, 🎏_Clay, hydro, ipClass, iZ, ksmodelτ, option, param, Ψ₁; 🎏_IsTopsoil=false, 🎏_RockFragment=false, RockFragment=[], IsTopsoil=[])

			return KsΨmodel = KSMODEL_OPTIONS(∑Psd, 🎏_Clay, 🎏_RockFragment, hydro, ipClass, iZ, ksmodelτ, option, param, Ψ₁; RockFragment=RockFragment)

		end  # function: KS_MODEL
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSMODEL_TRADITIONAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSMODEL_TRADITIONAL(Se, T1, T1Mac, T2, T2Mac, T3, T3Mac, θr, θs, θsMacMat, σ, σMac, Ψm, ΨmMac)

			Kunsat_Mat = T1 * ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2

			Kunsat_Mac = T1Mac * ((θs - θsMacMat) ^ T3Mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2Mac 

		return Kunsat_Mat + Kunsat_Mac
		end  # function: KS_MODEL
	# ------------------------------------------------------------------
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KsΨMODEL_NOINTEGRAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KsΨMODEL_NOINTEGRAL(T1, T1Mac, T2, T2Mac, T3, T3Mac, θr, θs, θsMacMat, σ, σMac, Ψ₁::Float64, Ψm, ΨmMac)

			# Matrix ====	
				θ_Mat = 0.5 * (θsMacMat - θr) * erfc((log( Ψ₁ / Ψm)) / (σ * √2.0)) + θr
				
				Se_Mat = (θ_Mat- θr) / (θs - θr)
				
				KsMat = T1 * cst.KunsatModel * π * ((θsMacMat - θr) * ((cst.Y / Ψm) ^ T2) * exp(((T2 * σ) ^ 2.0) / 2.0)) ^ T3
				
				Kunsat_Mat = KsMat * √Se_Mat * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0
			
			# Macropore ===
				θ_Mac = 0.5 * (θs - θsMacMat) * erfc((log(Ψ₁ / ΨmMac)) / (σMac * √2.0))
				
				Se_Mac = θ_Mac / (θs - θr)

				KsMac = T1Mac * cst.KunsatModel * π * ((θs - θsMacMat) * ((cst.Y / ΨmMac) ^ T2Mac) * exp(((T2Mac * σMac) ^ 2.0) / 2.0)) ^ T3Mac

				Kunsat_Mac = KsMac * √Se_Mac * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0
		return Kunsat_Mat + Kunsat_Mac
		end  # function: KsΨMODEL_NOINTEGRAL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KsΨMODEL_CLAY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KsΨMODEL_CLAY(hydro, iZ::Int64, optionₘ, T1, T1Mac, T2, T2Mac, T3, T3Mac, Tclay, θr, θs, θsMacMat, σ, σMac, Ψ₁::Float64, Ψm, ΨmMac)

			# Se ===
				Se = wrc.Ψ_2_Se(optionₘ, Ψ₁, iZ, hydro)

			# Matrix ====	
				KsMat = T1 * cst.KunsatModel * π * (((θsMacMat - θr) ^ (Tclay / T3)) * ((cst.Y / Ψm) ^ T2) * exp(((T2 * σ) ^ 2.0) / 2.0)) ^ T3

				Kunsat_Mat = KsMat * √Se * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0

			# Macropore ===
				KsMac = T1Mac * cst.KunsatModel * π * ((θs - θsMacMat) * ((cst.Y / ΨmMac) ^ T2Mac) * exp(((T2Mac * σMac) ^ 2.0) / 2.0)) ^ T3Mac

				Kunsat_Mac = KsMac * √Se * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0

	return  Kunsat_Mat + Kunsat_Mac
	end  # function: KsΨMODEL_CLAY
	# ------------------------------------------------------------------

	# ^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^__^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^__^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_^_

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TORTUOSITYMODELS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSMODEL_OPTIONS(∑Psd, 🎏_Clay::Bool, 🎏_RockFragment::Bool, hydro, ipClass, iZ::Int64, ksmodelτ, option, param, Ψ₁; RockFragment=[], θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψm=hydro.Ψm[iZ], σ=hydro.σ[iZ], θsMacMat=hydro.θsMacMat[iZ], ΨmMac=hydro.ΨmMac[iZ], σMac=hydro.σMac[iZ], τ₁ₐ=ksmodelτ.τ₁ₐ[ipClass],τclay₀=ksmodelτ.τclay₀[ipClass], τ₂ₐ=ksmodelτ.τ₂ₐ[ipClass], τclayₘₐₓ=ksmodelτ.τclayₘₐₓ[ipClass], τ₃ₐ=ksmodelτ.τ₃ₐ[ipClass], τclayΔθsr=ksmodelτ.τclayΔθsr[ipClass], τ₁ₐMac=ksmodelτ.τ₁ₐMac[ipClass],τclay₀Mac=ksmodelτ.τclay₀Mac[ipClass], τ₂ₐMac=ksmodelτ.τ₂ₐMac[ipClass], τclayₘₐₓMac=ksmodelτ.τclayₘₐₓMac[ipClass], τ₃ₐMac=ksmodelτ.τ₃ₐMac[ipClass], τclayΔθsrMac=ksmodelτ.τclayΔθsrMac)

			# Only correct if RF > Rf_StartIncrease
			if 🎏_RockFragment
				θr, θs, θsMacMat = ROCKCORRECTION!(hydro, iZ, RockFragment[iZ], θr, θs, θsMacMat)
			end #@isdefined RockFragment


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

					return Kunsat = KₛModel * kunsat.KUNSAT_θΨSe(option.hydro, Ψ₁, iZ, hydro) / hydro.Ks[iZ]
				
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
			return KsΨMODEL_NOINTEGRAL( T1, T1Mac, T2, T2Mac, T3, T3Mac, θr, θs, θsMacMat, σ, σMac, Ψ₁, Ψm, ΨmMac)
				
	
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
				
			else
				error("option.ksModel.KₛModel⍰ = $(option.ksModel.KₛModel⍰) is not yet implemented try <KsModel_Traditional>; <KsModel_Tσ>; <KsModel_New>; <KsModel_NewSimplified> ")
				
			end  # if: Model=="Model?"
		end  # function: TORTUOSITYMODELS 


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

				Clay = wrc.Ψ_2_Se(option.hydro, Ψ_Clay, iZ, hydro)
			end # 🎏_Clay
			
			X_Clay₁ =  τclay₀

			Clayₙ = max(Clay - X_Clay₁, 0.0) / (1.0 - X_Clay₁)

			ΔθsMacθr = hydro.θsMacMat[iZ] - hydro.θr[iZ]

			ΔθsMacθrₙ =  max(ΔθsMacθr - τclayΔθsr , 0.0) / (1.0 - τclayΔθsr)

			Tclay_Max =  1.0 + ΔθsMacθrₙ * (τclayₘₐₓ - 1.0) 

		return Tclay = Tclay_Max - (Tclay_Max - 1.0) * cos(Clayₙ * π * 0.5) 
		end				
	# ------------------------------------------------------------------

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ROCKCORRECTION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	""" The rock corection is already performed in θ(Ψ) and therefore Ks is already corected. Nevertheles, the model is wrong for RF > Rf_StartIncrease as the Ks starts to increase again"""
		function ROCKCORRECTION!(hydro, iZ, RockFragment, θr, θs, θsMacMat; Rf_StartIncrease=0.4, Rf_EndIncrease=0.9, θs_Amplify=1.1)

			Rf₁ = min(RockFragment, Rf_EndIncrease)

			if Rf₁ > Rf_StartIncrease
				# X values
					X = [Rf_StartIncrease, Rf_EndIncrease]
				
				# θs ----
					θs_NoRf = θs / (1.0 - Rf₁)
					Y_θs = [ (1.0 - Rf_StartIncrease) * θs_NoRf, θs_Amplify * θs_NoRf]
					Fit_θs = Polynomials.fit(X, Y_θs, 1)
					# θs = max(min(Fit_θs(Rf₁), hydro.θs_Max[iZ]), hydro.θs_Min[iZ])
					θs = Fit_θs(Rf₁)

				# θr ----
					θr_NoRf = θr / (1.0 - Rf₁)
					Y_θr = [(1.0 - Rf_StartIncrease) * θr_NoRf, θr_NoRf]
					Fit_θr = Polynomials.fit(X, Y_θr, 1)
					θr = max(min(Fit_θr(Rf₁), hydro.θr_Max[iZ]), hydro.θr_Min[iZ])

				# θsMacMat ----
					θsMacMat_NoRf =  θsMacMat / (1.0 - Rf₁)
					Y_θsMacMat = [min((1.0 - Rf_StartIncrease) * θsMacMat_NoRf, θs), 0.7 * (θs - θr) + θr]
					Fit_θsMacMat = Polynomials.fit(X, Y_θsMacMat, 1)	
					θsMacMat = min(Fit_θsMacMat(Rf₁), θs)					
			end
		return θr, θs, θsMacMat
		end  # function: ROCKCORRECTION
	# ------------------------------------------------------------------

end  # module θψ_2_KsψModel
