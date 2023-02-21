# =============================================================
#		module: kunsatModel jesus
# =============================================================
module θψ_2_KsψModel
	import ..cst, ..distribution, ..wrc
	import QuadGK
	import SpecialFunctions: erfc, erfcinv
	
	export KSΨMODEL_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KSΨMODEL_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KSΨMODEL_START(hydro, ipClass, iZ, ksmodelτ, KθModel, option, Ψ₁; Flag_IsTopsoil=false, Flag_RockFragment=false, IsTopsoil=[], RockFragment=[])

			if Flag_RockFragment
				RockFragment₁ = RockFragment[iZ]
			else
				RockFragment₁ = 0.0
			end #@isdefined RockFragment

			# if Flag_IsTopsoil
			# 	IsTopsoil₁ = Int64(IsTopsoil[iZ])
			# else
			# 	IsTopsoil₁ = 1	# Default value				
			# end  # if: @isdefined IsTopsoil

			KsΨmodel = TORTUOSITYMODELS(hydro, option, ipClass, iZ, ksmodelτ, Ψ₁; RockFragment=RockFragment₁, Smap_ImpermClass=[], KsImpClass_Dict=[])

		return KsΨmodel
		end  # function: KS_MODEL
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KsΨMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KsΨMODEL(hydro, iZ::Int64, optionₘ, T1, T1Mac, T2, T2Mac, T3, T3Mac, θr, θs, θsMacMat, σ, σMac, Ψ₁::Float64, Ψm, ΨmMac)

			# Se ===
				Se = wrc.kg.Ψ_2_SeDual(optionₘ, Ψ₁, iZ, hydro)

			# Matrix ====	
				Ks_Mat = T1 * cst.KunsatModel * π * ((θsMacMat - θr) * ((cst.Y / Ψm) ^ T2) * exp(((T2 * σ) ^ 2) / 2.0)) ^ T3

				Kunsat_Mat = Ks_Mat * √Se * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0

			# Macropore ===
				Ks_Mac = T1Mac * cst.KunsatModel * π * ((θs - θsMacMat) * ((cst.Y / ΨmMac) ^ T2Mac) * exp(((T2Mac * σMac) ^ 2) / 2.0)) ^ T3Mac

				Kunsat_Mac = Ks_Mac * √Se * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0
		return K_Ψ = Kunsat_Mat + Kunsat_Mac
		end  # function: KsΨMODEL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KsΨMODEL_OLD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		"""Pollacco, J.A.P., Webb, T., McNeill, S., Hu, W., Carrick, S., Hewitt, A., Lilburne, L., 2017. Saturated hydraulic conductivity model computed from bimodal water retention curves for a range of New Zealand soils. Hydrol. Earth Syst. Sci. 21, 2725–2737. https://doi.org/10.5194/hess-21-2725-2017"""
		function KsΨMODEL_OLD(hydro, iZ::Int64, optionₘ, T1, T1Mac, T2, T2Mac, T3, T3Mac, θr, θs, θsMacMat, σ, σMac, Ψ₁::Float64, Ψm, ΨmMac)

			# Ks Matrix ====	
				Ks_Mat = T1 * cst.KunsatModel * π * ((θsMacMat - θr) * ((cst.Y / Ψm) ^ T2) * exp(((T2 * σ) ^ 2) / 2.0)) ^ T3

			# Ks Macropore ===
				Ks_Mac = T1Mac * cst.KunsatModel * π * ((θs - θsMacMat) * ((cst.Y / ΨmMac) ^ T2Mac) * exp(((T2Mac * σMac) ^ 2) / 2.0)) ^ T3Mac
			
			Ks = Ks_Mat + Ks_Mac 
				
			# K(Ψ)
				Se = wrc.kg.Ψ_2_SeDual(optionₘ, Ψ₁, iZ, hydro)	

				 Kr_Mat = (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0

				 Kr_Mac = (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0

		return K_Ψ =  Ks * √Se * (Kr_Mat + Kr_Mac) 
		end  # function: KsΨMODEL_OLD
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TORTUOSITYMODELS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TORTUOSITYMODELS(hydro, option, ipClass, iZ::Int64,ksmodelτ, Ψ₁; RockFragment=0.0, θs=hydro.θs[iZ], θr=hydro.θr[iZ], Ψm=hydro.Ψm[iZ], σ=hydro.σ[iZ], θsMacMat=hydro.θsMacMat[iZ], ΨmMac=hydro.ΨmMac[iZ], σMac=hydro.σMac[iZ], τ₁ₐ=ksmodelτ.τ₁ₐ[ipClass],τ₁ᵦ=ksmodelτ.τ₁ᵦ[ipClass], τ₂ₐ=ksmodelτ.τ₂ₐ[ipClass], τ₂ᵦ=ksmodelτ.τ₂ᵦ[ipClass], τ₃ₐ=ksmodelτ.τ₃ₐ[ipClass], τ₃ᵦ=ksmodelτ.τ₃ᵦ[ipClass], τ₁ₐMac=ksmodelτ.τ₁ₐMac[ipClass],τ₁ᵦMac=ksmodelτ.τ₁ᵦMac[ipClass], τ₂ₐMac=ksmodelτ.τ₂ₐMac[ipClass], τ₂ᵦMac=ksmodelτ.τ₂ᵦMac[ipClass], τ₃ₐMac=ksmodelτ.τ₃ₐMac[ipClass], τ₃ᵦMac=ksmodelτ.τ₃ᵦMac, RockFragment_Treshold=0.4, Smap_ImpermClass=[], KsImpClass_Dict=[] )

			# Determine when Ks increases for increasing RockFragment	
				if RockFragment > RockFragment_Treshold
					θr, θs, θsMacMat = ROCKCORRECTION(RockFragment, RockFragment_Treshold, θr, θs, θsMacMat)
				end

				# MODEL 1 ====			
				if option.ksModel.KₛModel⍰=="KsΨmodel_1" # ===
					# Transformation matrix
						T1 = 10.0 ^ (τ₁ₐ / (τ₁ₐ - 1.0))
						T2_Min = 1.0; T2_Max = 3.0
						T2 = (T2_Min - T2_Max) * τ₂ₐ + T2_Max
						T3 = τ₃ₐ

					# Transformation macro
						T1Mac = 10.0 ^ (τ₁ₐMac / (τ₁ₐMac - 1.0))
						T2Mac = (T2_Min - T2_Max)  * τ₂ₐMac + T2_Max
						T3Mac = τ₃ₐMac									
					 return KsΨMODEL_OLD(hydro, iZ, option.hydro, T1, T1Mac, T2, T2Mac, T3, T3Mac, θr, θs, θsMacMat, σ, σMac, Ψ₁, Ψm, ΨmMac)

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
	#		FUNCTION :  τMODEL_σSilt
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_σSilt_σ(hydro, iZ, Pₘᵢₙ, Pₘₐₓ, σ; Amplitude=0.5, σSilt_η=0.538, Pσ=3, Distribution⍰="Normal", Normalise=true, Invert=false)
			ση = τMODEL_σ(hydro, iZ, Pₘₐₓ, Pₘᵢₙ, σ)

			τσ_Dist = τMODEL_σSilt(hydro, iZ, Pₘₐₓ, Pₘᵢₙ, ση; σSilt_η=σSilt_η, Pσ=Pσ, Distribution⍰="Normal", Normalise=Normalise, Invert=Invert)

			τσ = τMODEL_σ(hydro, iZ, Pₘₐₓ, Pₘᵢₙ, σ)

		return τ = min(τσ + Amplitude * (τσ_Dist / (σSilt_η + 1.0)) , 1.0) * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
		end
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  τMODEL_σSilt
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_σSilt(hydro, iZ, Pₘᵢₙ, Pₘₐₓ,  σ; σSilt_η=0.538, Pσ=3.0, Distribution⍰="Normal", Normalise=true, Invert=false)

			ση = τMODEL_σ(hydro, iZ, Pₘₐₓ, Pₘᵢₙ, σ)

			if  Distribution⍰== "Normal"
				σ_Dist = σSilt_η / Pσ

			elseif  Distribution⍰== "LogNormal"
				σ_Dist = log(σSilt_η) / Pσ

			else
				error("*** τMODEL_σSilt: $Distribution⍰ not implemented try <Normal> or  <LogNormal>  ***")
			end

			τσ_Dist = distribution.DISTRIBUTION(ση, σSilt_η, σ_Dist; Distribution⍰=Distribution⍰, Normalise=Normalise, Invert=Invert)[1]
		return τ = τσ_Dist  * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
		end
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : τMODEL_σ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_σ(hydro, iZ,  Pₘᵢₙ, Pₘₐₓ, σ; Inverse=false, τ₄=0.5)
			ση = σ_2_ση(hydro, iZ, σ)
			if Inverse
				return τ = (1.0 - ση) ^ τ₄  * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
			else
				return τ = ση ^ τ₄  * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
			end	
		end
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :τMODEL_σ2
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_σ2(hydro, iZ,  Pₘᵢₙ, Pₘₐₓ, σ; τ₄=0.5, σboundWater = 0.5)
			ση = σ_2_ση(hydro, iZ, σ)

		return τ = (Pₘₐₓ - Pₘᵢₙ) * min((ση / σboundWater) ^ τ₄, 1.0) + Pₘᵢₙ
		end
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : τMODEL_θsθr
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function τMODEL_θsθr(hydro, iZ,  Pₘᵢₙ, Pₘₐₓ, θs, θr, θsMacMat)
			θη = (θs - θr)
		return τ = (1.0 - θη)  * (Pₘₐₓ - Pₘᵢₙ) + Pₘᵢₙ
		end
	#..................................................................


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_ση
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function σ_2_ση(hydro, iZ, σ)
			return ση = (σ - hydro.σ_Min[iZ]) / (hydro.σ_Max[iZ] - hydro.σ_Min[iZ])
		end  # function: σ_2_ση
	# ------------------------------------------------------------------


end  # module θψ_2_KsψModel
# ............................................................