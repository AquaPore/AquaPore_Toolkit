	

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
			ση = σ_2_ση(hydro.σ[iZ])
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
		function σ_2_ση(Xσ)
			return ση = (Xσ - 0.75) / (3.5 - 0.75)
		end  # function: σ_2_ση
	# ------------------------------------------------------------------