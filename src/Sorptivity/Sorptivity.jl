# =============================================================
#		MODULE: sorptivity
# =============================================================
module sorptivity
	import ..wrc, ..kunsat
	import QuadGK
	import SpecialFunctions: erfc, erfcinv
	export SORPTIVITY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SORPTIVITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SORPTIVITY(θini::Float64, iZ::Int64, hydroInfilt₀, option, optionₘ; Rtol= 10^-8.0) #Rtol=10^-3.0

			# INITIALIZING

				Ψ_Sat = 0.0

				Se_Ini = wrc.θ_2_Se(θ₁=θini, θs=hydroInfilt₀.θs[iZ], θr=hydroInfilt₀.θr[iZ])

				SeIni_⬙ = (1.0 + Se_Ini) / 2.0

				θ⬙ = wrc.Se_2_θ(Se₁=SeIni_⬙, θs=hydroInfilt₀.θs[iZ], θr=hydroInfilt₀.θr[iZ])

				Ψ⬙ = wrc.θ_2_Ψ(optionₘ, θ⬙, iZ, hydroInfilt₀)

			# Sorptivity based on θ
				function SORPTIVITY_θ²(hydroInfilt₀, iZ, θ, θini, optionₘ)
					return DIFFUSIVITY_θ(θ, iZ, hydroInfilt₀, optionₘ) * (hydroInfilt₀.θs[iZ] + θ - 2.0 * θini)
				end # SORPTIVITY_θ² ~~~~~~~~~~~~~~~~~

				# Sorptivity_θ² = QuadGK.quadgk(θ -> SORPTIVITY_θ²(hydroInfilt₀, iZ, θ, θini, optionₘ), θini, θ⬙, rtol=Rtol)[1]

				# @show iZ
				# println("")
				# @show hydroInfilt₀.θs[iZ] 
				# @show hydroInfilt₀.θr[iZ]
				# @show hydroInfilt₀.Ks[iZ]
				# @show hydroInfilt₀.σ[iZ] 
				# @show hydroInfilt₀.Ψm[iZ]
				# @show hydroInfilt₀.θsMacMat[iZ]
				# @show hydroInfilt₀.ΨmMac[iZ] 
				@show hydroInfilt₀.σMac[iZ] 
		
				Sorptivity_θ² = QuadGK.quadgk(θ -> SORPTIVITY_θ²(hydroInfilt₀, iZ, θ, θini, optionₘ), θini, θ⬙)[1]

			# Sorptivity based on Ψ₁
				function SORPTIVITY_Ψ²(hydroInfilt₀, iZ, θini, Ψ₁)
					θ = wrc.Ψ_2_θ(optionₘ, Ψ₁, iZ, hydroInfilt₀)
				return kunsat.KUNSAT_θΨSe(optionₘ, Ψ₁, iZ, hydroInfilt₀) * (hydroInfilt₀.θs[iZ] + θ - 2.0 * θini)
				end # SORPTIVITY_Ψ² ~~~~~~~~~~~~~~~~~

				# Sorptivity_Ψ² = QuadGK.quadgk(Ψ₁ -> SORPTIVITY_Ψ²(hydroInfilt₀, iZ, θini, Ψ₁), Ψ_Sat, Ψ⬙, rtol=Rtol)[1]
				Sorptivity_Ψ² = QuadGK.quadgk(Ψ₁ -> SORPTIVITY_Ψ²(hydroInfilt₀, iZ, θini, Ψ₁), Ψ_Sat, Ψ⬙)[1]

		return √(max(Sorptivity_θ², eps()) + max(Sorptivity_Ψ², eps()))
		end # function SORPTIVITY


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#  	FUNCTION : DIFFUSIVITY
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function DIFFUSIVITY_θ(θ, iZ, hydroInfilt₀, optionₘ)

				Kunsat₀ = kunsat.KUNSAT_θΨSe(optionₘ, -1.0, iZ, hydroInfilt₀; θ₁=θ)

				# Ψ₁ = wrc.θ_2_Ψ(optionₘ, θ, iZ, hydroInfilt₀)

				Ψ₁ = wrc.θ_2_Ψ(optionₘ, θ, iZ, hydroInfilt₀)

			return - Kunsat₀ * wrc.∂Ψ∂θ(optionₘ, θ, iZ, hydroInfilt₀)
			# return - Kunsat₀ / wrc.∂θ∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroInfilt₀) #Faster but less stable
			end  # function: DIFFUSIVITY_θ ~~~~~~~~~~~~~~~~~


			# # DIFFUSIVITY_Se_η
			# function DIFFUSIVITY_Se_η(Se, iZ, hydroInfilt₀, optionₘ)
			# 	Kr = kunsat.Se_2_Kr(optionₘ, Se, iZ, hydroInfilt₀)

			# 	# Ψ₁ = wrc.Se_2_Ψ(optionₘ, Se, iZ, hydroInfilt₀)

			# return - Kr * wrc.∂Ψ∂Se(optionₘ, Se, iZ, hydroInfilt₀) / hydroInfilt₀.Ψm[iZ] # negative sign needed because Ψ₁ is set positive
			# end  # function: DIFFUSIVITY_Se_η ~~~~~~~~~~~~~~~~~


# 	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 	#		FUNCTION : FLUXXCONCENTRATION
# 	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 		function FLUXCONC(option, Se, Se_Ini)
# 			if option.infilt.SorptivityModel⍰== "Parlange" # <>=<>=<>=<>=<> soils with strong non-linear diffusivity behaviors
# 				return (1.0 + Se - 2.0 * Se_Ini) # to avoid numerical problems it is set analyticaly

# 			elseif option.infilt.SorptivityModel⍰== "Crank" # <>=<>=<>=<>=<> linear soil with constant diffusivity
# 				return (2.0 * (Se - Se_Ini)) / exp(-(erfcinv(max(Se - Se_Ini, eps()) / (1.0 - Se_Ini)))^2 )

# 			elseif option.infilt.SorptivityModel⍰== "Philip_Knight" # <>=<>=<>=<>=<> soils with Dirac δ-function diffusivity (i.e. Green and Ampt model)
# 				return (2.0 * (1.0 - Se_Ini)) # to avoid numerical problems it is set analyticaly 

# 			elseif option.infilt.SorptivityModel⍰== "Brutsaert" # <>=<>=<>=<>=<> soil with moderate non-linear diffusivity behaviors
# 				return (2.0 * √((Se - Se_Ini) * (1.0 - Se_Ini)))  # to avoid numerical problems it is set analyticaly
# 			end # option.infilt.SorptivityModel
# 		end  # function: FLUXXCONCENTRATION

end  # module: sorptivity
# ............................................................