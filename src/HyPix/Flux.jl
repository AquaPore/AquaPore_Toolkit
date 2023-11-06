module flux
	import ..cst
	import ..kunsat:KUNSAT_θΨSe 
	export Q!, K_AVER!

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Q
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Q!(optionHypix, discret, hydro, iZ::Int64, iT::Int64, Nz::Int64, paramHypix, Pkₐᵥₑᵣ::Vector{Float64}, Hpond, ΔPr_Soil, ΔSink, ΔT, θ, ψ_, ψ▲)
			if iZ == 1  # <>=<>=<>=<>=<>
				if optionHypix.TopBoundary⍰  == "Flux" 
					return ΔPr_Soil / ΔT[iT]

				elseif optionHypix.TopBoundary⍰ == "Ψ" 
					K_Aver₀, K_Aver = K_AVER!(optionHypix, paramHypix, Pkₐᵥₑᵣ, discret, hydro, iZ, Nz, ψ_, ψ▲)
					return K_Aver * (((ψ_ - paramHypix.Ψ_Top) / discret.ΔZ_⬓[1]) + paramHypix.Cosα)

				else
					error("Q! optionHypix.TopBoundary⍰ not found")
				end

			elseif 2 ≤ iZ ≤ Nz # <>=<>=<>=<>=<>
				K_Aver₀, K_Aver = K_AVER!(optionHypix, paramHypix, Pkₐᵥₑᵣ, discret, hydro, iZ, Nz, ψ_, ψ▲)
					return K_Aver * ( ((ψ_ - ψ▲) / discret.ΔZ_Aver[iZ]) + paramHypix.Cosα)

			else # iZ = Nz+1
				if optionHypix.BottomBoundary⍰ == "Free" # <>=<>=<>=<>=<>
					K_Aver₀, K_Aver = K_AVER!(optionHypix, paramHypix, Pkₐᵥₑᵣ, discret, hydro, iZ, Nz, ψ_, ψ▲)
						return K_Aver * paramHypix.Cosα

				elseif optionHypix.BottomBoundary⍰ == "Ψ" # <>=<>=<>=<>=<>
					K_Aver₀, K_Aver = K_AVER!(optionHypix, paramHypix, Pkₐᵥₑᵣ, discret, hydro, iZ, Nz, ψ_, ψ▲)
						return K_Aver * ( ((paramHypix.Ψ_Botom - ψ_) / discret.ΔZ_⬓[Nz]) + paramHypix.Cosα)
				else
					error("Q! optionHypix.BottomBoundary⍰ not found")
				end
			end # Case
		end  # function: Q!
	#-----------------------------------------------------------------



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : K_AVER!
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function K_AVER!(optionHypix, paramHypix, Pkₐᵥₑᵣ, discret, hydro, iZ::Int64, Nz::Int64, ψ_, ψ▲)
		
		K_Aver₀ = 1.0

		if iZ == 1 # <>=<>=<>=<>=<>
			if optionHypix.TopBoundary⍰  == "Flux" # <> = <> = <> = <> = <>
				K_Aver = 0.0::Float64

			elseif optionHypix.TopBoundary⍰  == "Ψ" # <> = <> = <> = <> = <>
				K_Aver = KUNSAT_θΨSe(optionHypix, ψ_, 1, hydro)

			else 	# <> = <> = <> = <> = <>
				error("K_AVER! optionHypix.TopBoundary⍰ not found")
			end	

		elseif 2 ≤ iZ ≤ Nz # <>=<>=<>=<>=<>
			K_Aver₀ = discret.ΔZ_W[iZ] * KUNSAT_θΨSe(optionHypix, ψ_, iZ, hydro) ^ Pkₐᵥₑᵣ[iZ] + (1.0 - discret.ΔZ_W[iZ]) * KUNSAT_θΨSe(optionHypix, ψ▲, iZ-1, hydro) ^ Pkₐᵥₑᵣ[iZ]

			K_Aver₀ = max(K_Aver₀, cst.Kθ_Min ^ Pkₐᵥₑᵣ[iZ])
			InvPkₐᵥₑᵣ = inv(Pkₐᵥₑᵣ[iZ])
			K_Aver = K_Aver₀ ^ InvPkₐᵥₑᵣ
			K_Aver₀ = K_Aver₀ ^ (InvPkₐᵥₑᵣ - 1.0)

		else # iZ = Nz+1
			K_Aver = KUNSAT_θΨSe(optionHypix, ψ_, Nz, hydro)
			K_Aver₀ = K_Aver
		end
	return K_Aver₀, K_Aver
	end  # function: K_AVER!
	#-------------------------------------------------------------------


	# =============================================================
	#		module: ∂Q∂ψ
	# 		only in use if ∂R∂Ψ_NumericalAuto" = false
	# =============================================================
	module ∂q∂Ψ
		import ..flux
		export ∂Q∂Ψ, ∂Q∂Ψ△, ∂Q▽∂Ψ, ∂Q▽∂Ψ▽
		import ...kunsat:KUNSAT_θΨSe 

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q∂Ψ(∂K∂Ψ::Vector{Float64}, discret, hydro, iT::Int64, iZ::Int64, K_Aver_Vect, K_Aver₀_Vect, Nz::Int64, optionHypix, paramHypix, Pkₐᵥₑᵣ, Ψ::Matrix{Float64})

				if iZ == 1 # <>=<>=<>=<>=<>
					if optionHypix.TopBoundary⍰ == "Flux" # =<>=<>=<>=<>=<>
						return 0.0::Float64

					elseif optionHypix.TopBoundary⍰ == "Ψ" # =<>=<>=<>=<>=<>		
						return ∂K∂Ψ[1] * ((Ψ[iT,1] - paramHypix.Ψ_Top) / discret.ΔZ_⬓[1] + paramHypix.Cosα) + K_Aver_Vect[1] / discret.ΔZ_⬓[1]
					
					else
						error("optionHypix.TopBoundary⍰ not found: ∂Q∂Ψ")
					end

				else # elseif 2 ≤ iZ ≤ Nz 	<>=<>=<>=<>=<>		
					# A1= discret.ΔZ_W[iZ] * ∂K∂Ψ[iZ] * ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + paramHypix.Cosα) + K_Aver / discret.ΔZ_Aver[iZ]

					return discret.ΔZ_W[iZ] * (KUNSAT_θΨSe(optionHypix, Ψ[iT,iZ], iZ, hydro) ^ (Pkₐᵥₑᵣ[iZ] - 1.0)) * K_Aver₀_Vect[iZ] * ∂K∂Ψ[iZ] * ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + paramHypix.Cosα) + K_Aver_Vect[iZ] / discret.ΔZ_Aver[iZ]	
				end # if iZ
			end  # function: ∂Q∂Ψ
		#-----------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q∂Ψ△
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q∂Ψ△(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, K_Aver_Vect, K_Aver₀_Vect, Nz::Int64, optionHypix, paramHypix, Pkₐᵥₑᵣ, Ψ)

				if iZ == 1 # <>=<>=<>=<>=<>
					return 0.0::Float64

				else # elseif 2 ≤ iZ ≤ Nz 	# <>=<>=<>=<>=<>
					# A1 =  (1.0 - discret.ΔZ_W[iZ]) * ∂K∂Ψ[iZ-1] *  ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + paramHypix.Cosα) - K_Aver / discret.ΔZ_Aver[iZ]

					return  (1.0 - discret.ΔZ_W[iZ]) * (KUNSAT_θΨSe(optionHypix, Ψ[iT,iZ-1], iZ, hydro) ^ (Pkₐᵥₑᵣ[iZ] - 1.0)) * K_Aver₀_Vect[iZ] * ∂K∂Ψ[iZ-1] * ((Ψ[iT,iZ] - Ψ[iT,iZ-1]) / discret.ΔZ_Aver[iZ] + paramHypix.Cosα) - K_Aver_Vect[iZ] / discret.ΔZ_Aver[iZ]
				end # if iZ
			end  # function: ∂Q∂Ψ△
		#-----------------------------------------------------------------

		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q▽∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q▽∂Ψ(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, K_Aver_Vect, K_Aver₀_Vect, Nz::Int64, optionHypix, paramHypix, Pkₐᵥₑᵣ, Ψ)
				
				if iZ ≤ Nz-1 	# <>=<>=<>=<>=<>
					# A1 =   (1.0 - discret.ΔZ_W[iZ+1]) * ∂K∂Ψ[iZ] * ((Ψ[iT,iZ+1] - Ψ[iT,iZ]) / discret.ΔZ_Aver[iZ+1] + paramHypix.Cosα) - K_Aver▽ / discret.ΔZ_Aver[iZ+1		

					return (1.0 - discret.ΔZ_W[iZ+1]) * ∂K∂Ψ[iZ] * (KUNSAT_θΨSe(optionHypix, Ψ[iT,iZ], iZ, hydro) ^ (Pkₐᵥₑᵣ[iZ+1] - 1.0)) * K_Aver₀_Vect[iZ+1] * ((Ψ[iT,iZ+1] - Ψ[iT,iZ]) / discret.ΔZ_Aver[iZ+1] + paramHypix.Cosα) - K_Aver_Vect[iZ+1] / discret.ΔZ_Aver[iZ+1]	
				
				else # iZ = Nz <>=<>=<>=<>=<>
					if optionHypix.BottomBoundary⍰ == "Free" # <>=<>=<>=<>=<>
						return ∂K∂Ψ[Nz] * paramHypix.Cosα
		
					elseif optionHypix.BottomBoundary⍰ == "Ψ" # <>=<>=<>=<>=<>
						return ∂K∂Ψ[Nz] * ((paramHypix.Ψ_Botom - Ψ[iT,Nz]) / discret.ΔZ_⬓[Nz] + paramHypix.Cosα) - K_Aver_Vect[Nz] /  discret.ΔZ_⬓[Nz]

					else
						error(" ∂Q▽∂Ψ optionHypix.BottomBoundary⍰ not found")
					end	
				end # if iZ
			end  # function: ∂Q▽∂Ψ
		#-----------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂Q▽∂Ψ▽
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂Q▽∂Ψ▽(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64, K_Aver_Vect, K_Aver₀_Vect, Nz::Int64, optionHypix, paramHypix, Pkₐᵥₑᵣ, Ψ)
				if iZ ≤ Nz-1 	# <>=<>=<>=<>=<>
					# A1 = discret.ΔZ_W[iZ+1] * ∂K∂Ψ[iZ+1] * ((Ψ[iT,iZ+1] - Ψ[iT,iZ]) / discret.ΔZ_Aver[iZ+1] + paramHypix.Cosα) + K_Aver▽ / discret.ΔZ_Aver[iZ+1]
					return discret.ΔZ_W[iZ+1] * (KUNSAT_θΨSe(optionHypix, Ψ[iT,iZ+1], iZ, hydro) ^ (Pkₐᵥₑᵣ[iZ+1] - 1.0)) * K_Aver₀_Vect[iZ+1] * ∂K∂Ψ[iZ+1] * ((Ψ[iT,iZ+1] - Ψ[iT,iZ]) / discret.ΔZ_Aver[iZ+1] + paramHypix.Cosα) + K_Aver_Vect[iZ+1] / discret.ΔZ_Aver[iZ+1]
					
				else # elseif iZ == Nz <>=<>=<>=<>=<>
					return 0.0::Float64
				end
			end  # function: ∂Q▽∂Ψ▽

	end  # module ∂q∂ψ
	# ............................................................

end # MODULE flux