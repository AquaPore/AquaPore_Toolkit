module kunsat
	import ..wrc
	export KUNSAT_θΨSe, Se_2_KUNSAT, θ_2_KUNSAT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KUNSAT_θΨSe
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KUNSAT_θΨSe(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			Ψ₁ = max(Ψ₁, 0.0)

			if  optionₘ.HydroModel⍰ == "Kosugi"
				return kunsat.kg.KUNSAT_θΨSe(Ψ₁=Ψ₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ], Option_KosugiModel_KΨ⍰=optionₘ.KosugiModel_KΨ⍰)

			elseif  optionₘ.HydroModel⍰ == "Vangenuchten" ||  optionₘ.HydroModel⍰ == "VangenuchtenJules"
				return kunsat.vg.KUNSAT_θΨSe(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "BrooksCorey"
				return kunsat.bc.KUNSAT_θΨSe(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "ClappHornberger"
				return kunsat.ch.KUNSAT_θΨSe(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for KUNSAT_θΨSe is not yet available")
			end
		end # function KUNSAT_θΨSe
	#-------------------------------------------------------------------

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θ_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function θ_2_KUNSAT(optionₘ, θ₁, iZ::Int64, hydroParam)
			
			if  optionₘ.HydroModel⍰ == "Kosugi"
				return kunsat.kg.KUNSAT_θΨSe(θ₁=θ₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ],  Option_KosugiModel_KΨ⍰=optionₘ.KosugiModel_KΨ⍰)
			else
				error("$( optionₘ.HydroModel⍰) model for θ_2_KUNSAT is not yet available")
			end
		end # function θ_2_KUNSAT
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function Se_2_KUNSAT(optionₘ, Se₁, iZ::Int64, hydroParam)
			Se = max(min(Se₁, 1.0), 0.0)

			if  optionₘ.HydroModel⍰ == "Kosugi"

				return kunsat.kg.KUNSAT_θΨSe(Se₁=Se₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ], Option_KosugiModel_KΨ⍰=optionₘ.KosugiModel_KΨ⍰)
				
			elseif  optionₘ.HydroModel⍰ == "Vangenuchten"
				return kunsat.vg.Se_2_KUNSAT(optionₘ, Se₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "BrooksCorey"
				return kunsat.bc.Se_2_KUNSAT(optionₘ, Se₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "ClappHornberger"
				return kunsat.ch.Se_2_KUNSAT(optionₘ, Se₁, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for Se_2_KUNSAT is not yet available")
			end
		end # function Se_2_KUNSAT
	#-------------------------------------------------------------------
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_Kr
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function Se_2_Kr(optionₘ, Se₁, iZ::Int64, hydroParam)
			Se₁ = max(min(Se₁, 1.0), 0.0)

			if  optionₘ.HydroModel⍰ == "Kosugi"
				return  kunsat.kg.KUNSAT_θΨSe(Se₁=Se₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ], Option_KosugiModel_KΨ⍰=optionₘ.KosugiModel_KΨ⍰) / hydroParam.Ks[iZ]
			else
				error("$( optionₘ.HydroModel⍰) model for Se_2_Kr is not yet available")
			end
		end # function Se_2_Kr
	#-------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂K∂ΨMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function ∂K∂ΨMODEL(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			Ψ₁ = max(Ψ₁, 0.0)

			if  optionₘ.HydroModel⍰ == "Kosugi"
				return kunsat.kg.∂K∂ΨMODEL(Ψ₁=Ψ₁, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

			elseif  optionₘ.HydroModel⍰ == "Vangenuchten"
				return kunsat.vg.∂K∂ΨMODEL(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "BrooksCorey"
				return kunsat.bc.∂K∂ΨMODEL(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			elseif  optionₘ.HydroModel⍰ == "ClappHornberger"
				return kunsat.ch.∂K∂ΨMODEL(optionₘ, Ψ₁, iZ::Int64, hydroParam)
			else
				error("$( optionₘ.HydroModel⍰) model for ∂K∂ΨMODEL is not yet available")
			end
		end # function ∂K∂ΨMODEL
	#-------------------------------------------------------------------


	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# =============================================================
	#		MODULE KOSUGI
	# =============================================================
	module kg
		import..wrc
		import ...cst, ...hydroRelation
		import ForwardDiff, QuadGK
		import SpecialFunctions: erfc, erfcinv
		export KUNSAT_θΨSe, ∂K∂ΨMODEL

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KUNSAT_θΨSe
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function KUNSAT_θΨSe(;Ψ₁=-1.0, θ₁=-1.0, Se₁ =-1.0, θs, θsMacMat, θr, Ψm, σ, ΨmMac, σMac, Ks, Option_KosugiModel_KΨ⍰="Traditional")

				if Ψ₁==-1.0 && θ₁==-1.0 && Se₁==-1.0
					error("KUNSAT_θΨSe function: Cannot 3 of them: Ψ₁==-1.0 && θ₁=-1.0 && Se₁=-1.0 ")
				end

				# For more flexibility on the different inputs
					if θ₁ ≠ -1.0 && Se₁ == -1.0
						Se₁ = wrc.θ_2_Se(θ₁=θ₁, θs=θs, θr=θr)

					elseif  Se₁ == -1.0 && Ψ₁ ≠ -1.0				
						Se₁ = wrc.kg.Ψ_2_Se(Ψ₁=Ψ₁, θs=θs, θsMacMat=θsMacMat, θr=θr, Ψm=Ψm, σ=σ, ΨmMac=ΨmMac, σMac=σMac)

					else
						error("KUNSAT_θΨSe function: needs: Se₁ or θ data")

					end 

					if Ψ₁ == -1.0 && θ₁ ≠ -1.0
						Ψ₁ = wrc.kg.θ_2_Ψ(θ₁=θ₁, θs=θs, θsMacMat=θsMacMat, θr=θr, Ψm=Ψm, σ=σ, ΨmMac=ΨmMac, σMac=σMac) 

					elseif Ψ₁ == -1.0 && Se₁ ≠ -1.0
						Ψ₁ = wrc.kg.Se_2_Ψ(Se₁=Se₁, θs=θs, θsMacMat=θsMacMat, θr=θr, Ψm=Ψm, σ=σ, ΨmMac=ΨmMac, σMac=σMac) 

					end

				if  Option_KosugiModel_KΨ⍰ == "Traditional" # =====

					KsMat = Ks * min(max((θsMacMat - θr) / (θs - θr), 0.0), 1.0)			
					Kunsat_Mat =  KsMat * √Se₁ * (0.5 * erfc(((log(Ψ₁/ Ψm)) / σ + σ) / √2.0)) ^ 2.0

					KsMac = max(Ks - KsMat, 0.0)
					Kunsat_Mac =  KsMac * √Se₁ * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0

					return Kunsat_Mat + Kunsat_Mac
				

				elseif Option_KosugiModel_KΨ⍰ == "Mualem" # =====

					Kunsat_Mat = 0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)
			
					Kunsat_Mac = 0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)

					W_Mat = (θsMacMat - θr) * exp((σ^2.0) / 2.0) / Ψm 
		
					W_Mac = max(θs - θsMacMat, 0.0) * exp((σMac^2.0) / 2.0) / ΨmMac

					return Ks * √Se₁ * ((W_Mat * Kunsat_Mat + W_Mac * Kunsat_Mac) / (W_Mat + W_Mac)) ^ 2.0
				

				elseif Option_KosugiModel_KΨ⍰ == "ΨmacMat" # =====
				
					ΨmacMat = hydroRelation.FUNC_θsMacMatη_2_ΨmacMat(;θs, θsMacMat, θr)

					KsMat = Ks * min(max((θsMacMat - θr) / (θs - θr), 0.0), 1.0)

					Se_Mat = 0.5 * erfc((log( max(Ψ₁ - ΨmacMat, 0.0) / Ψm)) / (σ * √2.0))

					Kunsat_Mat =  KsMat * √Se_Mat * (0.5 * erfc(((log( max(Ψ₁- ΨmacMat, 0.0)/ Ψm)) / σ + σ) / √2.0)) ^ 2.0

					KsMac = max(Ks - KsMat, 0.0)
					Kunsat_Mac =  KsMac * √Se₁ * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0

					return Kunsat_Mat + Kunsat_Mac

				elseif Option_KosugiModel_KΨ⍰ == "Quirijn" # =====

					Pq = 0.5

					Kunsat_Mat = (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0
			
					Kunsat_Mac = (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0

					W_Mat₀ = ((θsMacMat - θr) * exp((σ^2.0) / 2.0) / Ψm ) ^ Pq
		
					W_Mac₀ = (max(θs - θsMacMat, 0.0) * exp((σMac^2.0) / 2.0) / ΨmMac) ^ Pq

					W_Mat = W_Mat₀ / (W_Mac₀ + W_Mat₀)
					
					W_Mac = W_Mac₀ / (W_Mac₀ + W_Mat₀)

					return Ks * √Se₁ * (W_Mat * Kunsat_Mat + W_Mac * Kunsat_Mac)
			
				else
					error("option.hydro.Option_KosugiModel_KΨ⍰ = $Option_KosugiModel_KΨ⍰ not yet available pls modify ?_Option.toml")
				end

			end # function KUNSAT_θΨSe
		#-------------------------------------------------------------------


		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# #		FUNCTION : θ_2_KUNSAT
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	function θ_2_KUNSAT(;θ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, σMac, Ks)

		# 		Se₁ = wrc.θ_2_Se(θ₁=θ₁, θs=θs, θr=θr)

		# 		Ψ₁ = wrc.kg.θ_2_Ψ(θ₁=θ₁, θs=θs, θsMacMat=θsMacMat, θr=θr, Ψm=Ψm, σ=σ, ΨmMac=ΨmMac, σMac=σMac) 
	
		# 		KsMat = Ks * min(max((θsMacMat - θr) / (θs - θr), 0.0), 1.0)			
		# 		Kunsat_Mat =  KsMat * √Se₁ * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0

		# 		KsMac = max(Ks - KsMat, 0.0)
		# 		Kunsat_Mac =  KsMac * √Se₁ * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0

		# 	return Kunsat_Mat + Kunsat_Mac		
		# 	end # function θ_2_KUNSAT
		# #-------------------------------------------------------------------

		
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# #		FUNCTION : Se_2_KUNSAT
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	function Se_2_KUNSAT(;Se₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, σMac, Ks)

		# 		Ψ₁ = wrc.kg.Se_2_Ψ(Se₁=Se₁, θs=θs, θsMacMat=θsMacMat, θr=θr, Ψm=Ψm, σ=σ, ΨmMac=ΨmMac, σMac=σMac) 
	
		# 		KsMat = Ks * min(max((θsMacMat - θr) / (θs - θr), 0.0), 1.0)			
		# 		Kunsat_Mat = KsMat * √Se₁ * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0

		# 		KsMac = max(Ks - KsMat, 0.0)
		# 		Kunsat_Mac =  KsMac * √Se₁ * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0

		# 	return Kunsat_Mat + Kunsat_Mac			
		# 	end # function Se_2_KUNSAT
		# #-------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂K∂ΨMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂ΨMODEL(;Ψ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, σMac, Ks)

				if Ψ₁ > eps(100.0)
					SeMat = (θsMacMat - θr) / (θs - θr)
					SeMac =(θs - θsMacMat) / (θs - θr)

					KsMat = Ks * SeMat 
					KsMac = Ks * SeMac

					P₂ = 0.7071067811865475 # 1.0 / √(2.0) 
					Pπ = 0.5641895835477563 # 1.0 / √(π)
					Pc = 0.3989422804014327 # 1.0 / √(2.0*π) 
					P8 = 0.125 # 1.0 / 8.0

					Ψmσ = Ψm * σ
					ΨmMacσMac = ΨmMac * σMac
					ΨmΨ₁ = Ψm / Ψ₁
					ΨmMacΨ₁ = ΨmMac / Ψ₁
					Ψ₁Ψm = Ψ₁ / Ψm
					Ψ₁ΨmMac = Ψ₁ / ΨmMac
					F_LOG = (P₂ * log(Ψ₁Ψm)) / σ
					Erfc1_Mat = erfc(P₂ * σ + F_LOG)
					ΨΨm = erfc(log(Ψ₁Ψm) / (√2.0 *σ))
					ΨΨmσ = exp((-(log(Ψ₁Ψm) / (√2.0 *σ))*log(Ψ₁Ψm)) / (√2.0 *σ))
					F_LOG_MAC = (P₂ * log(Ψ₁ΨmMac)) / σMac
					Erfc1_Mac = erfc(P₂ * σMac + F_LOG_MAC)
					ΨΨmMac = erfc(log(Ψ₁ΨmMac) / (√2.0 *σMac))
					ΨΨmσ_Mac = exp((-(log(Ψ₁ΨmMac) / (√2.0 *σMac))*log(Ψ₁ΨmMac)) / (√2.0 *σMac))
					Sqrt = sqrt(0.5*SeMat*ΨΨm+0.5*SeMac*ΨΨmMac)
				
					∂Kunsat_Mat∂Ψ = (-Pc * KsMat*(ΨmΨ₁)*Erfc1_Mat*Sqrt*exp(((-F_LOG) - P₂*σ) * (P₂*σ + F_LOG))) / (Ψmσ)+ P8*KsMat*(-Pπ)*((SeMac*(ΨmMacΨ₁)*ΨΨmσ_Mac) / (√2.0 *ΨmMacσMac) + (SeMat*(ΨmΨ₁)*ΨΨmσ) / (√2.0 *Ψmσ))*(Erfc1_Mat^2)*(Sqrt^-1)
					
					if θs - θsMacMat > cst.ΔθsθsMacMat 
						∂Kunsat_Mac∂Ψ =  (-Pc * KsMac*(ΨmMacΨ₁)*Erfc1_Mac*Sqrt*exp((-F_LOG_MAC - P₂*σMac)*(P₂ * σMac + F_LOG_MAC))) / (ΨmMacσMac) + P8 *(-Pπ)*KsMac*((SeMac*(ΨmMacΨ₁)*ΨΨmσ_Mac) / (√2.0 *ΨmMacσMac) + (SeMat*(ΨmΨ₁)*ΨΨmσ) / (√2.0 *Ψmσ))*(Erfc1_Mac^2)*(Sqrt^-1.0)
					else
						∂Kunsat_Mac∂Ψ = 0.0::Float64
					end
					
					return ∂Kunsat_Mac∂Ψ + ∂Kunsat_Mat∂Ψ
				else
					return 0.0::Float64
				end 
			end
		#--------------------------------------------------------------------------------


		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# #		FUNCTION : ∂K∂SE
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	function ∂K∂SE(optionₘ, Ψ₁, Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])
				
		# 		if Se < eps(100.0)
		# 			Se += eps(10.0)
		# 		end

		# 		P1 = 1.0 / sqrt(2.0)
		# 		P2 = 0.125

		# 		KsMat = Ks * (θsMacMat - θr) / (θs - θr)

		# 		∂Kunsat_Mat∂θ = KsMat * √Se * exp(-(σ / √2.0 + erfcinv( 2.0*Se )) ^ 2 + erfcinv(2.0*Se )^2.0)*erfc(σ/ √2.0 + erfcinv(2.0*Se)) + P2*KsMat*erfc(σ / √2.0 + erfcinv(2.0*Se)) ^2 / √Se

		# 		if θs - θsMacMat > cst.ΔθsθsMacMat
		# 			KsMac = Ks * (θs - θsMacMat) / (θs - θr)

		# 			∂Kunsat_Mac∂θ = KsMac*sqrt(Se)*exp(-(σMac/ √2.0 + erfcinv(2.0*Se)) ^ 2 + erfcinv(2.0*Se) ^2.0)*erfc(σMac / √2.0 + erfcinv(2.0*Se)) + P2*KsMac*erfc(P1*σMac + erfcinv(2.0*Se)) ^ 2 / √Se
		# 		else
		# 			∂Kunsat_Mac∂θ = 0.0
		# 		end

		# 	return ∂Kunsat_Mat∂θ + ∂Kunsat_Mac∂θ
		# 	end #  ∂K∂θ
		# #--------------------------------------------------------------------------------

		

		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# #		FUNCTION : ∂K∂ΨMODEL numerical
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	function ∂K∂Ψ_NUMERICAL(optionₘ, Ψ₁, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])
				
		# 		ψ =fill(0.0::Float64, 1) 

		# 		∂K∂Ψ_Numerical(ψ::Vector) = KUNSAT_θΨSe(optionₘ, abs(ψ[1]), iZ, hydroParam)
				
		# 		ψ[1] = Ψ₁
				
		# 		Func_∂K∂Ψ_Numerical = ψ -> ForwardDiff.gradient(∂K∂Ψ_Numerical, ψ)			
		# 		∂K∂Ψ = Func_∂K∂Ψ_Numerical(ψ)[1]
				
		# 		if isnan(∂K∂Ψ)
		# 			error("∂K∂Ψ = NaN")
		# 			∂K∂Ψ = 0.0
		# 		end
		# 	return ∂K∂Ψ 
		# 	end # function: ∂K∂Ψ

			
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# #		FUNCTION : ∂K∂Ψ analitical
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	function ∂K∂Ψ_Jesus3(optionₘ, Ψ₁, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

		# 		Se = wrc.Ψ_2_Se(optionₘ, Ψ₁, iZ, hydroParam)

		# 		KsMat = Ks * (θsMacMat - θr) / (θs - θr)
		# 		KsMac = Ks * (θs - θsMacMat)  / (θs - θr)

		# 		Ψ₁Ψm = log(Ψ₁ / Ψm)
		# 		Ψ₁ΨmMac = log(Ψ₁ / ΨmMac)

		# 		A = Ψ₁Ψm / (σ * √2.0)
		# 		B = Ψ₁ΨmMac / (σMac * √2.0)

		# 	  ∂Kunsat_Mat∂Ψ = -(KsMat / (√(2.0*π)*Ψ₁*σ)) * (erfc(A+σ/√2.0)) * √Se * (exp(-(A+σ/√2.0)^2.0)) - (KsMat / 8*(√(2.0*π)*Ψ₁)) * ((erfc(A+σ/√2.0))^2.0) / √Se * ((θsMacMat-θr)/(θs-θr) * exp(-(A)^2.0)/σ + (θs-θsMacMat)/(θs-θr) * exp(-(B)^2.0)/σMac) 
				
		# 	  if isnan(∂Kunsat_Mat∂Ψ)
		# 			∂Kunsat_Mat∂Ψ = 0.0
		# 		end

		# 		∂Kunsat_Mac∂Ψ = -(KsMac / (√(2.0*π)*Ψ₁*σMac)) * (erfc(B+σMac/√2.0)) * √Se * (exp(-(B+σMac/√2.0)^2.0)) - (KsMat / 8*(√(2.0*π)*Ψ₁)) * ((erfc(B+σMac/√2.0))^2.0) / √Se * ((θsMacMat-θr)/(θs-θr) * exp(-(A)^2.0)/σ + (θs-θsMacMat)/(θs-θr) * exp(-(B)^2.0)/σMac)

		# 		if isnan(∂Kunsat_Mac∂Ψ)
		# 			∂Kunsat_Mac∂Ψ = 0.0
		# 		end

		# 	return ∂Kunsat_Mat∂Ψ + ∂Kunsat_Mac∂Ψ
		# 	end # function: ∂K∂ΨMODEL analitical
		# #-----------------------------------------------------------------

				# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# #		FUNCTION : KUNSAT_θΨSe
			# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# 	function Ψ_2_KUNSAT2(;Ψ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, σMac, Ks)

			# 		Se = wrc.kg.Ψ_2_Se(Ψ₁=Ψ₁, θs=θs, θsMacMat=θsMacMat, θr=θr, Ψm=Ψm, σ=σ, ΨmMac=ΨmMac, σMac=σMac)

			# 		ΨmacMat = exp(log(ΨmMac) + 3.0 * σMac)

			# 		SeMacMat = max((θsMacMat - θr) / (θs - θr), 0.0)

			# 		# Kr_ΨmacMat =  1.0 + √SeMacMat * (0.5 * erfc(((log(ΨmacMat / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0

			# 		Kr_ΨmacMat =  1.0 + √SeMacMat 
					
			# 		KsMac = Ks * min(1.0 / Kr_ΨmacMat, 1.0)
			# 		KsMat = Ks - KsMac
					
			# 		Kunsat_Mac =  KsMac * √Se * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0

			# 		Kunsat_Mat =  KsMat * √Se * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0

			# 	return Kunsat_Mat + Kunsat_Mac
			# 	end # function KUNSAT_θΨSe
			# #-------------------------------------------------------------------
	end # module kg



	# =============================================================
	#		MODULE VAN GENUCHTEN
	# =============================================================
	module vg
		import ..wrc
		export KUNSAT_θΨSe, Se_2_KUNSAT, ∂K∂ΨMODEL

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : KUNSAT_θΨSe
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  KUNSAT_θΨSe(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Ks=hydroParam.Ks[iZ], Km=hydroParam.Km[iZ], L=0.5)
				M = 1.0 - Km / N

				Se = wrc.Ψ_2_Se(optionₘ, Ψ₁, iZ, hydroParam)
				return Kunsat = Ks * (Se^L) * ( 1.0 - (1.0 - Se ^ (1.0 / M) ) ^ M ) ^ 2.0
			end #function KUNSAT_θΨSe


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Se_2_KUNSAT(optionₘ, Se, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Ks=hydroParam.Ks[iZ], Km=hydroParam.Km[iZ], L=0.5)
				M = 1.0 - Km / N
			return Kunsat = Ks * (Se.^L) * ( 1. - (1. - Se.^ (1.0 / M) ) .^ M ) .^ 2.0
			end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : ∂K∂ΨMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂ΨMODEL(optionₘ, Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Ks=hydroParam.Ks[iZ], Km=hydroParam.Km[iZ], L=0.5)
				M = 1.0 - Km/N
		
				∂K∂Ψ = Ks * (L * (-M * (N * (1 / Ψvg) * (Ψ₁ / Ψvg) ^ (N - 1)) * (1.0 + (Ψ₁ / Ψvg) ^ N) ^ (-M - 1)) * ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (L - 1)) * (1.0 - (1.0 - ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ M) ^ 2.0 + Ks *
				((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ L * (2.0 * -(M * -((1.0 / M) * (-M * (N * (1 / Ψvg) * (Ψ₁ / Ψvg) ^ (N - 1)) * (1.0 + (Ψ₁ / Ψvg) ^ N) ^ (-M - 1)) * ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M - 1)) * (1.0 - ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ (M - 1)) * (1.0 - (1.0 - ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ M))

			return ∂K∂Ψ
			end # function ∂K∂Ψ

	end #module vg ...............................................


	# =============================================================
	#		MODULE BROOKS AND COOREY
	# =============================================================
	module bc
		import ..wrc
		export KUNSAT_θΨSe, Se_2_KUNSAT, ∂K∂ΨMODEL

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : KUNSAT_θΨSe
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  KUNSAT_θΨSe(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λbc - 2.0

				if Ψ₁ > Ψbc
					return Kunsat = Ks * (Ψ₁ / Ψbc) ^ M 
				else 
					return Kunsat = Ks
				end
			end #function KUNSAT_θΨSe


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Se_2_KUNSAT(optionₘ, Se, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λbc - 2.0

			return Ks * Se .^ M 
			end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : ∂K∂ΨMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂ΨMODEL(optionₘ, Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λbc - 2.0
		
			return Ks * M / Ψbc * (Ψ₁/Ψbc) ^ (M-1.0) 
			end # function ∂K∂Ψ


	end #module bc ...............................................


	# =============================================================
	#		MODULE CLAPP AND HORNBERGER
	# =============================================================
	module ch
		import ..wrc
		export KUNSAT_θΨSe, Se_2_KUNSAT, ∂K∂ΨMODEL

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : KUNSAT_θΨSe
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  KUNSAT_θΨSe(optionₘ, Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λch - 2.0

				if Ψ₁ > Ψch
					return Ks * (Ψ₁ / Ψch) ^ M 
				else 
					return Ks
				end
			end #function KUNSAT_θΨSe


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Se_2_KUNSAT(optionₘ, Se, iZ, hydroParam, θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψbc[iZ], λch=hydroParam.λch[iZ], Ks=hydroParam.Ks[iZ])
					
				M = -3.0 * λch - 2.0
			return Ks * Se .^ M 
			end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : ∂K∂ΨMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂ΨMODEL(optionₘ, Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λch - 2.0
		
			return Ks * M / Ψch * (Ψ₁/Ψch) ^ (M-1.0) 
			end # function ∂K∂ΨMODEL

	end #module ch ...............................................


end # module kunsat 
# ...........................................................................