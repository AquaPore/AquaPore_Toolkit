# =============================================================
#		MODULE WRC
# =============================================================
module wrc
	export ∂Se∂Ψ, ∂θ∂Ψ, ∂Ψ∂Se, ∂Ψ∂θ, Se_2_θ, θ_2_Se, θ_2_Ψ, Ψ_2_Se, Ψ_2_θ

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θ_2_Se
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θ_2_Se(;θ₁, θs, θr)
			Se = (θ₁ - θr) / (θs - θr)
		return max( min(Se, 1.0), 0.0)
		end # function θ_2_Se
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_θ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Se_2_θ(;Se₁, θs, θr)
			θ₁ = Se₁ * (θs - θr) + θr
		return max(min(θ₁, θs), θr)
		end # function Se_2_θ
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_θ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_θ(optionₘ, Ψ₁, iZ, hydroParam)
			Ψ₁ = max(Ψ₁, 0.0)

			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.Ψ_2_θ(Ψ₁=Ψ₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], ΨmacMat=hydroParam.ΨmacMat[iZ], σMac=hydroParam.σMac[iZ], KosugiModel_θΨ⍰=optionₘ.KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac=optionₘ.ΨmacMat_2_σMac_ΨmMac)

			elseif optionₘ.HydroModel⍰ == "Vangenuchten"
				return wrc.vg.Ψ_2_θ(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "VangenuchtenJules"
				return wrc.vgJules.Ψ_2_θ(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "BrooksCorey"
				return wrc.bc.Ψ_2_θ(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "ClappHornberger"
				return wrc.ch.Ψ_2_θ(Ψ₁, iZ::Int64, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for Ψ_2_θ is not yet available")
			end
		end # function Ψ_2_θ
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_Se
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_Se(optionₘ, Ψ₁, iZ, hydroParam)
			Ψ₁ = max(Ψ₁, 0.0)

			if optionₘ.HydroModel⍰ == "Kosugi"

				return wrc.kg.Ψ_2_Se(Ψ₁=Ψ₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], ΨmacMat=hydroParam.ΨmacMat[iZ], σMac=hydroParam.σMac[iZ], KosugiModel_θΨ⍰=optionₘ.KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac=optionₘ.ΨmacMat_2_σMac_ΨmMac)

			elseif optionₘ.HydroModel⍰ == "Vangenuchten"
				return wrc.vg.Ψ_2_Se(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "VangenuchtenJules"
				return wrc.vgJules.Ψ_2_Se(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "BrooksCorey"
				return wrc.bc.Ψ_2_Se(Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "ClappHornberger"
				return wrc.ch.Ψ_2_Se(Ψ₁, iZ::Int64, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for Ψ_2_Se is not yet available")
			end # function Ψ_2_Se
		end
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θ_2_Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θ_2_Ψ(optionₘ, θ₁, iZ, hydroParam)
			if optionₘ.HydroModel⍰ == "Kosugi"

				return wrc.kg.θ_2_Ψ(θ₁=θ₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], ΨmacMat=hydroParam.ΨmacMat[iZ], σMac=hydroParam.σMac[iZ], KosugiModel_θΨ⍰=optionₘ.KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac=optionₘ.ΨmacMat_2_σMac_ΨmMac)

			elseif optionₘ.HydroModel⍰ == "Vangenuchten"
				return wrc.vg.θ_2_Ψ(θ₁, iZ, hydroParam)
			elseif optionₘ.HydroModel⍰ == "BrooksCorey"
				return wrc.bc.θ_2_Ψ(θ₁, iZ, hydroParam)
			elseif optionₘ.HydroModel⍰ == "ClappHornberger"
				return wrc.ch.θ_2_Ψ(θ₁, iZ, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for θ_2_Ψ is not yet available")
			end # function θ_2_Ψ
		end
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Se_2_Ψ(optionₘ, Se₁, iZ, hydroParam)
			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.Se_2_Ψ(Se₁=Se₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], ΨmacMat=hydroParam.ΨmacMat[iZ], σMac=hydroParam.σMac[iZ], KosugiModel_θΨ⍰=optionₘ.KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac=optionₘ.ΨmacMat_2_σMac_ΨmMac)
			else
				error("$(optionₘ.HydroModel⍰) model for Se_2_Ψ is not yet available")
			end # function Se_2_Ψ
		end
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂θ∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂θ∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			Ψ₁ = max(Ψ₁, 0.0)

			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.∂θ∂Ψ(Ψ₁=Ψ₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], ΨmacMat=hydroParam.ΨmacMat[iZ], σMac=hydroParam.σMac[iZ], KosugiModel_θΨ⍰=optionₘ.KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac=optionₘ.ΨmacMat_2_σMac_ΨmMac)

			elseif optionₘ.HydroModel⍰ == "Vangenuchten"
				return wrc.vg.∂θ∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "BrooksCorey"
				return wrc.bc.∂θ∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "ClappHornberger"
				return wrc.ch.∂θ∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for ∂θ∂Ψ is not yet available")
			end
		end # function ∂θ∂Ψ
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂Ψ∂θ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂Ψ∂θ(optionₘ, θ₁, iZ::Int64, hydroParam)

			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.∂Ψ∂θ(θ₁=θ₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], ΨmacMat=hydroParam.ΨmacMat[iZ], σMac=hydroParam.σMac[iZ], KosugiModel_θΨ⍰=optionₘ.KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac=optionₘ.ΨmacMat_2_σMac_ΨmMac)
			else
				error("$(optionₘ.HydroModel⍰) model for ∂Ψ∂θ is not yet available")
			end
		end # function ∂Ψ∂θ
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂Se∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂Se∂Ψ(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			Ψ₁ = max(Ψ₁, 0.0::Float64)

			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.∂Se∂Ψ(Ψ₁=Ψ₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], ΨmacMat=hydroParam.ΨmacMat[iZ], σMac=hydroParam.σMac[iZ], KosugiModel_θΨ⍰=optionₘ.KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac=optionₘ.ΨmacMat_2_σMac_ΨmMac)
			else
				error("$(optionₘ.HydroModel⍰) model for ∂Se∂Ψ is not yet available")
			end
		end # function ∂Se∂Ψ
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂Se∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂Ψ∂Se(optionₘ, Se₁, iZ::Int64, hydroParam)

			if optionₘ.HydroModel⍰ == "Kosugi"
				return wrc.kg.∂Ψ∂Se(Se₁=Se₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], ΨmacMat=hydroParam.ΨmacMat[iZ], σMac=hydroParam.σMac[iZ], KosugiModel_θΨ⍰=optionₘ.KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac=optionₘ.ΨmacMat_2_σMac_ΨmMac)
			else
				error("$(optionₘ.HydroModel⍰) model for ∂Ψ∂Se is not yet available")
			end
		end # function ∂Se∂Ψ
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : GREEN-AMPT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function GREEN_AMPT(optionₘ, iZ::Int64, hydroParam)
			if optionₘ.HydroModel⍰ == "BrooksCorey"
				return  wrc.bc.GREEN_AMPT(optionₘ, iZ::Int64, hydroParam)
			elseif optionₘ.HydroModel⍰ == "ClappHornberger"
				return wrc.ch.GREEN_AMPT(optionₘ, iZ::Int64, hydroParam)
			else
				error("$(optionₘ.HydroModel⍰) model for GREEN_AMPT is not yet available")
			end
		end # function: GREEN_AMPT


	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# =============================================================
	#		MODULE KOSUGI
	# =============================================================
	module kg
		import SpecialFunctions: erfc, erfcinv
		import ForwardDiff, Optim
		import ...cst, ...hydroRelation, ..wrc
		export ∂Se∂Ψ, ∂θ∂Ψ, ∂Ψ∂Se, ∂Ψ∂θ, Se_2_Ψ, θ_2_Ψ, Ψ_2_Se, Ψ_2_θ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_θ(;Ψ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac=100.0, ΨmacMat, σMac=1.1, KosugiModel_θΨ⍰="Traditional", ΨmacMat_2_σMac_ΨmMac=true, Pσ_Mac=2.0)

				# Physically constraining the hydraulic parameters
				if ΨmacMat_2_σMac_ΨmMac == true
					ΨmMac = hydroRelation.FUNC_ΨmacMat_2_ΨmMac(;ΨmacMat)
					σMac  = hydroRelation.FUNC_ΨmacMat_2_σMac(;ΨmacMat)
				end

				if KosugiModel_θΨ⍰ == "Traditional" || KosugiModel_θΨ⍰ == "Mualem"
					θ_Mat = 0.5 * (θsMacMat - θr) * erfc((log( Ψ₁ / Ψm)) / (σ * √2.0)) + θr

					θ_Mac = 0.5 * max(θs - θsMacMat, 0.0) * erfc((log(Ψ₁ / ΨmMac)) / (σMac * √2.0))
					return θ_Mac + θ_Mat


				elseif KosugiModel_θΨ⍰ == "Traditional_Constrained"
						θ_Mat = 0.5 * (θsMacMat - θr) * erfc((log( Ψ₁ / Ψm)) / (σ * √2.0)) + θr

						θ_Mac = 0.5 * max(θs - θsMacMat, 0.0) * erfc((Pσ_Mac / √2.0 ) * (log(Ψ₁) / log(√ΨmacMat) - 1.0))
						return θ_Mac + θ_Mat


				elseif KosugiModel_θΨ⍰ == "ΨmacMat"
					if Ψ₁ ≤ ΨmacMat
						return θ_Mac =  max(θs - θsMacMat, 0.0) * 0.5 * (erfc((log(Ψ₁ / ΨmMac)) / (σMac * √2.0)) - (min(Ψ₁ / ΨmacMat, 1.0)) * erfc((log(ΨmacMat / ΨmMac)) / (σMac * √2.0))) + θsMacMat
					else
						return θ_Mat = 0.5 * (θsMacMat - θr) * erfc((log( max(Ψ₁ - ΨmacMat, eps(100.0)) / Ψm)) / (σ * √2.0)) + θr
					end

				else
					error("function Ψ_2_θ has no $KosugiModel_θΨ⍰")
				end
			end # function Ψ_2_θ
		#-----------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_Se(;Ψ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰="Traditional", ΨmacMat_2_σMac_ΨmMac=true)

				θ₁ = wrc.kg.Ψ_2_θ(;Ψ₁, θs, θsMacMat,θr, Ψm, σ=σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac)

			return wrc.θ_2_Se(;θ₁, θs, θr)
			end # function Ψ_2_Se
		#-----------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θ_2_Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θ_2_Ψ(;θ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰="Traditional", ΨmacMat_2_σMac_ΨmMac=true)

				function OF(Ψ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, σMac)
					θmod = wrc.kg.Ψ_2_θ(;Ψ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac)
				return (θ₁ - θmod) ^ 4.0
				end # OF

				if θs - θsMacMat < 0.0001
					Se = max( min((θ₁ - θr) / (θs - θr), 1.0) , 0.0)
					Ψ₀ = Ψm * exp(erfcinv(2.0 * Se) * σ * √2.0)
				else
					Optimization = Optim.optimize(Ψ₁ -> OF(10.0 ^ Ψ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, σMac), log10(0.0001), log10(1.0E6), Optim.GoldenSection())
					Ψ₀ = 10.0 ^ Optim.minimizer(Optimization)[1]
				end
			return max(Ψ₀, 0.0)
			end # θ_2_Ψ
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Se_2_Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Se_2_Ψ(;Se₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰="Traditional", ΨmacMat_2_σMac_ΨmMac=true)
				θ₁ = wrc.Se_2_θ(;Se₁, θs, θr)
			return  θ_2_Ψ(;θ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac)
			end # Se_2_Ψ
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂θ∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂Ψ(;Ψ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰="Traditional", ΨmacMat_2_σMac_ΨmMac=true)

				# Physically constraining the hydraulic parameters
				if ΨmacMat_2_σMac_ΨmMac == true
					ΨmMac = hydroRelation.FUNC_ΨmacMat_2_ΨmMac(;ΨmacMat)
					σMac  = hydroRelation.FUNC_ΨmacMat_2_σMac(;ΨmacMat)
				end

				if Ψ₁ > eps(100.0)
					∂θ∂Ψ_Mat = (θsMacMat - θr) * exp( -((log(Ψ₁ / Ψm)) ^ 2.0) / (2.0 * σ ^ 2.0))  / (Ψ₁ * σ * √(π * 2.0))

					∂θ∂Ψ_Mac = (θs - θsMacMat) * exp( -((log(Ψ₁ / ΨmMac)) ^ 2.0) / (2.0 * σMac ^ 2.0)) / (Ψ₁ * σMac * √(π * 2.0))

					return ∂θ∂Ψ_Mat + ∂θ∂Ψ_Mac
				else
					return 0.0::Float64
				end
			end # function ∂θ∂Ψ
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂θ∂Ψ_NORM
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂Ψ_NORM(;Ψ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰="Traditional", ΨmacMat_2_σMac_ΨmMac=true)

				# Physically constraining the hydraulic parameters
				if ΨmacMat_2_σMac_ΨmMac == true
					ΨmMac = hydroRelation.FUNC_ΨmacMat_2_ΨmMac(;ΨmacMat)
					σMac  = hydroRelation.FUNC_ΨmacMat_2_σMac(;ΨmacMat)
				end

				# ΨmacMat = hydroRelation.FUNC_θsMacMatη_2_ΨmacMat(;θs, θsMacMat, θr)

				Ψmod_Mat = exp(log(Ψm) - σ^2)

				Ψmod_Mac = exp(log(ΨmMac) - σMac^2)

				if Ψ₁ > eps(100.0)
					∂θ∂Ψ_Mat(Ψ₁) = (θsMacMat - θr) * exp( -((log(max(Ψ₁-ΨmacMat , 0.0) / Ψm)) ^ 2.0) / (2.0 * σ^2.0)) / (Ψ₁ * σ * √(π * 2.0))

					∂θ∂Ψ_Mat_Mod = (θsMacMat - θr) * exp( -((log(Ψmod_Mat / Ψm)) ^ 2.0) / (2.0 * σ^2.0)) / (Ψmod_Mat * σ * √(π * 2.0))

					∂θ∂Ψ_Mac(Ψ₁) = (θs - θsMacMat) * exp( -((log(Ψ₁ / ΨmMac)) ^ 2.0) / (2.0 * σMac^2.0)) / (Ψ₁ * σMac * √(π * 2.0))

					∂θ∂Ψ_Mac_Mod = (θs - θsMacMat) * exp( -((log(Ψmod_Mac / ΨmMac)) ^ 2.0) / (2.0 * σMac^2.0)) / (Ψmod_Mac * σMac * √(π * 2.0))

					return ∂θ∂Ψ_Mat(Ψ₁) / ∂θ∂Ψ_Mat(Ψmod_Mat) + ∂θ∂Ψ_Mac(Ψ₁) / ( ∂θ∂Ψ_Mac(Ψmod_Mac) + eps())
				else
					return 0.0
				end # function ∂θ∂Ψ_NORM
			end
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂θ∂R
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂R(;R₁, θs, θsMacMat, θr, Rm, σ, RmMac, σMac, KosugiModel_θΨ⍰="Traditional", ΨmacMat_2_σMac_ΨmMac=true)

				if R₁ > eps()
					∂θ∂R_Mat = ((θsMacMat - θr) * exp( -((log(R₁ / Rm)) ^ 2.0) / (2.0 * σ ^ 2.0))) / (R₁ * σ * √(π * 2.0))

					∂θ∂R_Mac = ((θs - θsMacMat) * exp( -((log(R₁ / RmMac)) ^ 2.0) / (2.0 * σMac ^ 2.0))) / (R₁ * σMac * √(π * 2.0))

					return ∂θ∂R_Mat + ∂θ∂R_Mac
				else
					return 0.0::Float64
				end
			end # function ∂θ∂R
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂θ∂R_NORM
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂R_NORM(;R₁, θs, θsMacMat, θr, Rm, σ, RmMac, σMac, KosugiModel_θΨ⍰="Traditional", ΨmacMat_2_σMac_ΨmMac=true)


				# ΨmacMat = hydroRelation.FUNC_θsMacMatη_2_ΨmacMat(;θs, θsMacMat, θr)

				RMacMat = cst.Y / ΨmacMat

				if R₁ <  eps(100.0)
					R₁ +=  eps()
				end

				∂θ∂R_Mat(R₁) = -(θsMacMat - θr) * exp( -((log(max(R₁ - RMacMat,0.0)/ Rm)) ^ 2.0) / (2.0 * σ ^ 2.0)) / (R₁ * σ * √(π * 2.0))

				∂θ∂R_Mac(R₁) = - max(θs - θsMacMat, 0.0) * exp( -((log(R₁ / RmMac)) ^ 2.0) / (2.0 * σMac ^ 2.0)) / (R₁ * σMac * √(π * 2.0))

				Rmod_Mat = exp(log(Rm) - σ^2)

				Rmod_Mac = exp(log(RmMac) - σMac^2)

			return ∂θ∂R_Mat(R₁) / ∂θ∂R_Mat(Rmod_Mat) + ∂θ∂R_Mac(R₁) / ∂θ∂R_Mac(Rmod_Mac)
			end # function ∂θ∂R_NORM
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Ψ∂θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Ψ∂θ(;θ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰="Traditional", ΨmacMat_2_σMac_ΨmMac=true)

				θ₁ = max(min(θs - eps(1000.0), θ₁), θr + eps(1000.0))

				θ₂ = fill(0.0::Float64, 1)

				∂Ψ∂θ_Numerical(θ₂::Vector) = θ_2_Ψ(θ₁=abs(θ₂[1]), θs=θs, θsMacMat=θsMacMat, θr=θr, Ψm=Ψm, σ=σ, ΨmMac=ΨmMac, σMac=σMac, ΨmacMat=ΨmacMat, KosugiModel_θΨ⍰=KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac=ΨmacMat_2_σMac_ΨmMac)

				θ₂[1] = θ₁
				Func_∂Ψ∂θ_Numerical = θ₂ -> ForwardDiff.gradient(∂Ψ∂θ_Numerical, θ₂)

			return ∂Ψ∂θ₀ = Func_∂Ψ∂θ_Numerical(θ₂)[1]

				# if isnan(∂Ψ∂θ₀)
				# 	println( "======== wrc.kg.∂Ψ∂θ = NaN ======== ")
				# 	return 1.0
				# else
				# 	return ∂Ψ∂θ₀
				# end
			end # function ∂Ψ∂θ
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Ψ∂Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Ψ∂Se(;Se₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰="Traditional", ΨmacMat_2_σMac_ΨmMac=true)

				Se₁ = max(min(Se₁, 1.0 - eps(1000.0)), eps(1000.0))

				Se₂ = fill(0.0::Float64, 1)
				# ∂Ψ∂Se_Numerical(Se₂::Vector) = Se_2_Ψ(optionₘ, abs(Se₂[1]), iZ, hydroParam)
				∂Ψ∂Se_Numerical(Se₂::Vector) = Se_2_Ψ(Se₁=abs(Se₂[1]), θs=θs, θsMacMat=θsMacMat, θr=θr, Ψm=Ψm, σ=σ, ΨmMac=ΨmMac, σMac=σMac, KosugiModel_θΨ⍰=KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac=ΨmacMat_2_σMac_ΨmMac)

				Se₂[1] = Se₁
				Func_∂Ψ∂Se_Numerical = Se₂ -> ForwardDiff.gradient(∂Ψ∂Se_Numerical, Se₂)

				# ∂Ψ∂Se₀ = Func_∂Ψ∂Se_Numerical(Se₂)[1]
			return Func_∂Ψ∂Se_Numerical(Se₂)[1]
			end # function ∂Se∂Ψ
		#-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Se∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Se∂Ψ(;Ψ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰="Traditional", ΨmacMat_2_σMac_ΨmMac=true)

				Ψ₁ = max(Ψ₁, eps(1000.0))
				ψ = fill(0.0::Float64, 1)

				# ∂Ψ∂Se_Numerical(ψ::Vector) = Ψ_2_Se(optionₘ, abs(ψ[1]), iZ, hydroParam)
				 ∂Ψ∂Se_Numerical(ψ::Vector) = Ψ_2_Se(Ψ₁= abs(ψ[1]), θs=θs, θsMacMat=θsMacMat, θr=θr, Ψm=Ψm, σ=σ, ΨmMac=ΨmMac, σMac=σMac, KosugiModel_θΨ⍰=KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac=ΨmacMat_2_σMac_ΨmMac)
				ψ[1] = Ψ₁
				Func_∂Ψ∂Se_Numerical = ψ -> ForwardDiff.gradient(∂Ψ∂Se_Numerical, ψ)

			return Func_∂Ψ∂Se_Numerical(ψ)[1]
			end # function ∂Ψ∂Se


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂θ∂Ψ Mode
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# function ∂θ∂Ψ_Mode(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])

			# 	Ψm_Mode = exp(log(Ψm)-σ^2)

			# 	∂θ∂Ψ_Mode = ∂θ∂Ψ(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=Ψm_Mode, σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ])

			# return Ψm_Mode, ∂θ∂Ψ_Mode
			# end # function ∂θ∂Ψ_Mode
		#-------------------------------------------------------------------

	end # module kg # ...............................................



	# ===============================================================================================
	#		MODULE VAN GENUCHTEN
	# ===============================================================================================
	module vg
		import ..wrc
		export Ψ_2_θ, Ψ_2_Se, θ_2_Ψ, ∂θ∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : Ψ_2_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ]) # van Genuchten WRC
				M = 1.0 - Km / N
				Se = (1.0 + (Ψ₁ / Ψvg) ^ N ) ^ (-M)
			return θ₂ = wrc.Se_2_θ(Se₁=Se, θs=θs, θr=θr)
			end # function Ψ_2_θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_Se(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ]) # van Genuchten WRC
				M = 1.0 - Km / N
			return Se = (1.0 + (Ψ₁ / Ψvg) ^ N ) ^ (-M)
			end # function Ψ_2_Se


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : θ_2_Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θ_2_Ψ(θ₂, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ]) # van Genuchten WRC
				θ₂ = max(min(θ₂, θs), θr)

				Se = wrc.θ_2_Se(θ₁=θ₂, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ])

				M = 1. - Km / N
			return Ψ₁ = Ψvg * exp(log(exp(log(Se) / -M) - 1.0) / N)
			end # Se_2_Ψ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : ∂θ∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂Ψ(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ])
				M = 1.0 - Km/N # van Genuchten
			return M * (θs-θr) / ( Ψvg*(1-M)) * ((Ψ₁/Ψvg) .^ (N * M)) .* (1.0 + (Ψ₁/Ψvg) .^ N) .^ (-M-1.)
			end # function ∂θ∂Ψ

	end # module vg # ..................................................................................................................



	# ======================================================================================
	#		MODULE BROOKS AND COREY
	# ======================================================================================
	module bc
		import ..wrc
		export Ψ_2_θ, ∂θ∂Ψ, Ψ_2_Se, θ_2_Ψ, GREEN_AMPT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : Ψ_2_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ]) # Brooks & Corey WRC
				if Ψ₁ > Ψbc
					return θ₂ = θr + (θs - θr) * (Ψ₁ / Ψbc) ^ - λbc
				else
					return θ₂ = θs
				end
			end # function Ψ_2_θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_Se(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ])

				θ₂ = Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ])

				return Se = wrc.θ_2_Se(θ₁=θ₂, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ])
			end # function Ψ_2_Se


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : θ_2_Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θ_2_Ψ(θ₂, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ])

				Se = wrc.θ_2_Se(θ₁=θ₂, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ])
				return Ψ₁ = Ψbc * (Se ^ -λbc)
			end # θ_2_Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : ∂θ∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂Ψ(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ])

				return (-λbc / Ψbc) * (θs-θr) * (Ψ₁/Ψbc) .^ (-λbc - 1.)
			end # function ∂θ∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : GREEN_AMPT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			""" Parameters derived by Brooks and Corey 1964
			Brooks RH, Corey AT. 1964. Hydraulic properties of porous media. Hydrology Papers 3, Colorado State University, Fort Collins, 27 p. 3"""
			function GREEN_AMPT(optionₘ, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], A=2.0, B=3.0, C=1.0)
				return Ψga = Ψbc * (A + B *λbc) / (C + B * λbc)
			end # GREEN_AMPT

		end # module brooksCorey # ...............................................


	# ======================================================================================
	#		MODULE CLAPP AND HORNBERGER
	# ======================================================================================
	module ch
		import ..wrc
		export Ψ_2_θ, ∂θ∂Ψ, Ψ_2_Se, θ_2_Ψ, GREEN_AMPT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch: Ψ_2_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ]) # Clapp & Hornberger WRC
				if Ψ₁ > Ψch
					return θ₂ = θs * (Ψ₁ / Ψch) ^ - λch
				else
					return θ₂ = θs
				end
			end # function Ψ_2_θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_Se(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ])

				θ₂ = Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ])

			return wrc.θ_2_Se(θ₁=θ₂, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ])
			end # function Ψ_2_Se

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : θ_2_Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θ_2_Ψ(θ₂, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ])

				Se = wrc.θ_2_Se(θ₁=θ₂, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ])
			return Ψch * (Se ^ -λch)
			end # θ_2_Ψ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : ∂θ∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂θ∂Ψ(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ])

			return (-λch / Ψch) * θs * (Ψ₁/Ψch) .^ (-λch - 1.)
			end # function ∂θ∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : GREEN_AMPT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			""" Parameters derived by Brooks and Corey 1964
			Brooks RH, Corey AT. 1964. Hydraulic properties of porous media. Hydrology Papers 3, Colorado State University, Fort Collins, 27 p. 3"""
			function GREEN_AMPT(optionₘ, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ], A=2.0, B=3.0, C=1.0)
			return Ψga = Ψch * (A + B *λch) / (C + B * λch)
			end # GREEN_AMPT

	end # module CLAPP AND HORNBERGER # ...............................................



	# ===============================================================================================
	#		MODULE VAN GENUCHTEN JULES
	# ===============================================================================================
	module vgJules
		import ..wrc
		export Ψ_2_θ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		vg FUNCTION : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_θ(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ]) # van Genuchten WRC
				M = 1.0 - Km / N
				if Ψ₁ > Ψvg
					Se = ( (Ψ₁ / Ψvg) ^ N ) ^ (-M)
				else
					Se =1.0
				end
			return θ₂ = wrc.Se_2_θ(Se₁=Se, θs=θs, θr=θr)
			end # function Ψ_2_θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_Se
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_Se(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Km=hydroParam.Km[iZ]) # van Genuchten WRC
				M = 1.0 - Km / N
				if Ψ₁ > Ψvg
					Se = (1+ (Ψ₁ / Ψvg) ^ N ) ^ (-M)
				else
					Se = 1.0
				end
			end # function Ψ_2_See

	end # module vgJules ...............................................

end # module wrc # ...............................................
