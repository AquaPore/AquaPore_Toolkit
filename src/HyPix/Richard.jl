# =============================================================
#		MODULE: residual
# =============================================================
module richard
	import ..timeStep, ..flux, ..ponding, ..residual
	import ..wrc: Ψ_2_θDual, ∂θ∂Ψ, θ_2_ΨDual
	import ..kunsat
	using LinearAlgebra

	export RICHARD_ITERATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD_ITERATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge::Bool, Hpond::Vector{Float64}, hydro, iCount_ReRun::Int64, iNonConverge::Int64, iT::Int64, IterCount::Int64, K_Aver_Vect::Vector{Float64}, K_Aver₀_Vect::Vector{Float64}, Nz::Int64, optionHypix, paramHypix, Pkₐᵥₑᵣ::Vector{Float64}, Q, Residual, Sorptivity::Float64, ΔLnΨmax::Vector{Float64}, ΔPrThroughfall::Vector{Float64}, ΔRunoff::Vector{Float64}, ΔSink::Matrix{Float64}, ΔT, θ::Matrix{Float64}, Ψ_Max::Float64, Ψ_Min::Vector{Float64}, Ψ::Matrix{Float64}, Ψbest::Vector{Float64})
						
			# INITIALIZING
			 @simd for iZ = 1:Nz
					Ψ[iT,iZ] = Ψ[iT-1,iZ]
				end # for iZ = 1:Nz
	
			# ITTERATION
			Residual_Max_Best = Inf
			iTer = 0::Int64
			ΔPr_Soil = 0.0::Float64
			while iTer ≤ paramHypix.N_Iter - 1
            iTer      += 1
            IterCount += 1 # Counting the iterations

				# RESIDUAL MAX BEST: Deriving the Residual max because may be Ψ[iT-1,iZ] is the best solution
				∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Hpond, Q, Residual, ΔPr_Soil, θ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, Hpond, hydro, iT, K_Aver_Vect, K_Aver₀_Vect, Nz::Int64, optionHypix, paramHypix, Pkₐᵥₑᵣ, Q, Residual, Sorptivity, ΔPrThroughfall, ΔRunoff, ΔSink, ΔT, θ, Ψ, Ψ_Max)

				# Computing Residual_Max_Best at the beginning before iteration
				if iTer == 1
					Residual_Max_Best = CONVERGENCECRITERIA(discret, iT, Nz, Residual, ΔT)
				end

				# The minimum Ψ depends if we are close to saturation
					Ψ_Min = ΨMIN(iT, Nz, paramHypix, Ψ, Ψ_Min)

					Ψ = SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT, Nz, optionHypix, paramHypix, Residual, ΔLnΨmax, θ, Ψ, Ψ_Min, Ψ_Max)

				# Averaging the Residuals, depending on method
					Residual_Max = CONVERGENCECRITERIA(discret, iT, Nz,  Residual, ΔT)

				# Determine if iteration made improvement
					if Residual_Max < Residual_Max_Best	
						@simd for iZ=1:Nz
							Ψbest[iZ] = Ψ[iT,iZ]
						end
						Residual_Max_Best = copy(Residual_Max)
					end # Residual_Max < Residual_Max_Best 	

				# Did we achieve the goals
				if Residual_Max ≤ paramHypix.WaterBalanceResidual_Max
					break # Move out the loop
				end  # if: Residual
			end # while: iTer ======================

			# Making sure we get the best if convergence fails
			# No convergence
			if iTer == paramHypix.N_Iter
				Flag_NoConverge = true

				# Put the best values
				@simd for iZ=1:Nz
					Ψ[iT,iZ] = Ψbest[iZ]
				end
			else
				Flag_NoConverge = false
			end #  iTer == paramHypix.N_Iter

			# UPDATE Θ
			@simd for iZ=1:Nz
				θ[iT,iZ] = Ψ_2_θDual(optionHypix, Ψ[iT,iZ], iZ, hydro)
			end

			# Determine if the simulation is going to rerun with a different time step
			Flag_ReRun, iCount_ReRun, iNonConverge, ΔT = RERUN_HYPIX(discret, Flag_NoConverge, Hpond, hydro, iCount_ReRun, iNonConverge, iT, Nz, optionHypix, paramHypix, Pkₐᵥₑᵣ, Q, ΔLnΨmax, ΔPrThroughfall, ΔPr_Soil, ΔSink, ΔT, θ, Ψ)

		return Flag_NoConverge, Flag_ReRun, Hpond, iCount_ReRun, iNonConverge, iTer, IterCount, Q, ΔRunoff, ΔT, θ, Ψ
		end  # function: RICHARD_SOLVING
	#----------------------------------------------------------]-------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge::Bool, Hpond, hydro, iT::Int64, K_Aver_Vect, K_Aver₀_Vect, Nz::Int64, optionHypix, paramHypix, Pkₐᵥₑᵣ, Q, Residual, Sorptivity, ΔPrThroughfall, ΔRunoff, ΔSink, ΔT, θ, Ψ, Ψ_Max)
		
			Hpond, ΔPr_Soil, ΔRunoff = ponding.PONDING_RUNOFF_SORPTIVITY(discret, Hpond, hydro, iT, optionHypix, paramHypix, Sorptivity, ΔPrThroughfall, ΔRunoff, ΔSink, ΔT, θ, Ψ)
		
		#----------------------------------------------------------------
			# ∂R∂Ψ2 = fill(0.0, Nz)
			# ∂R∂Ψ▽2 = fill(0.0, Nz)
			# ∂R∂Ψ△2 =  fill(0.0, Nz)

			for iZ=1:Nz
				Q, Residual, θ = residual.RESIDUAL(discret, hydro, iT, iZ, Nz, optionHypix, paramHypix, Pkₐᵥₑᵣ, Q, Residual, Hpond, ΔPr_Soil, ΔSink, ΔT, θ, Ψ)

				∂K∂Ψ[iZ] = kunsat.∂K∂ΨMODEL(optionHypix, Ψ[iT,iZ], iZ, hydro)
				
				K_Aver₀_Vect[iZ], K_Aver_Vect[iZ] = flux.K_AVER!(optionHypix, paramHypix, Pkₐᵥₑᵣ, discret, hydro, iZ, Nz, Ψ[iT,iZ], Ψ[iT,max(iZ-1,1)])
			end
			K_Aver₀_Vect[Nz+1], K_Aver_Vect[Nz+1] = flux.K_AVER!(optionHypix, paramHypix, Pkₐᵥₑᵣ, discret, hydro, Nz+1, Nz, Ψ[iT,Nz], Ψ[iT,Nz])

			for iZ=1:Nz
				if !(optionHypix.∂R∂Ψ_NumericalAuto)
					∂R∂Ψ[iZ], ∂R∂Ψ△[iZ], ∂R∂Ψ▽[iZ] = residual.∂RESIDUAL∂Ψ(∂K∂Ψ, discret, hydro, iT, iZ, Nz, optionHypix, paramHypix, Pkₐᵥₑᵣ, K_Aver_Vect, K_Aver₀_Vect, ΔT, θ, Ψ)
				else
					∂R∂Ψ[iZ] = residual.∂R∂Ψ_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, Nz, optionHypix, paramHypix, Pkₐᵥₑᵣ, Hpond, ΔPr_Soil, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,Nz)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,Nz)], Ψ_Max)

					∂R∂Ψ▽[iZ]  = residual.∂R∂Ψ▽_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, Nz, optionHypix, paramHypix, Pkₐᵥₑᵣ, Hpond, ΔPr_Soil, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,Nz)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,Nz)], Ψ_Max)
		
					∂R∂Ψ△[iZ]  = residual.∂R∂Ψ△_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, Nz, optionHypix, paramHypix, Pkₐᵥₑᵣ, Hpond, ΔPr_Soil, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,Nz)], Ψ[iT,iZ], Ψ[iT, min(iZ+1,Nz)], Ψ_Max)
					
				end # if optionHypix.∂R∂Ψ_NumericalAuto"
			end #for iZ= 1:Nz

			# # FOR TESTING...
				# println("One:=================")
				# println("∂R∂Ψ_Deriv=" , ∂R∂Ψ[1:Nz],"\n") 
				# println("∂R∂Ψ_Num=" , ∂R∂Ψ2[1:Nz],"\n")
				# println("∂R∂Ψ_Deriv=" , ∂R∂Ψ[1:Nz] .- ∂R∂Ψ2[1:Nz],"\n") # No good at cell N

				# println("Two: =================")
				# println("∂R∂Ψ▽_Num=" , ∂R∂Ψ▽[1:Nz],"\n")
				# println("∂R∂Ψ▽_Der=" , ∂R∂Ψ▽2[1:Nz],"\n") # No good
				# println("∂R∂Ψ▽_Num=" , ∂R∂Ψ▽[1:Nz] .- ∂R∂Ψ▽2[1:Nz],"\n")

				# println("Tree: =================")
				# println("∂R∂Ψ△_Num=" , ∂R∂Ψ△[1:Nz],"\n") # Good
				# println("∂R∂Ψ△_Der=" , ∂R∂Ψ△2[1:Nz],"\n")
				# println("∂R∂Ψ△_Der=" , ∂R∂Ψ△[1:Nz] .- ∂R∂Ψ△2[1:Nz],"\n")

		return ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Hpond, Q, Residual, ΔPr_Soil, θ
		end # function RICHARD
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΨMIN
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ΨMIN(iT::Int64, Nz::Int64, paramHypix, Ψ::Matrix{Float64}, Ψ_Min::Vector{Float64})

			Flag_Ψsmall = false
			for iZ=1:Nz
				if Ψ[iT-1,iZ] < paramHypix.opt.ΨmacMat / 2.0 # mm
					Flag_Ψsmall = true
					break
				end
			end

			@simd for iZ=1:Nz
				if Flag_Ψsmall
					Ψ_Min[iZ] = paramHypix.Ψ_MinMin
				else
					Ψ_Min[iZ] = 0.0::Float64
				end
			end

		# 	@simd for iZ=1:Nz
		# 		if Ψ[iT-1,iZ] < paramHypix.opt.ΨmacMat / 2.0 # mm
		# 			Ψ_Min[iZ] = paramHypix.Ψ_MinMin
		# 		else
		# 			Ψ_Min[iZ] = 0.0::Float64
		# 		end
		# 	end
		return Ψ_Min
		end  # function: ΨMIN
	# ------------------------------------------------------------------



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERGENCECRITERIA
	#     Averaging the Residuals, depending on method
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CONVERGENCECRITERIA(discret, iT::Int64, Nz::Int64, Residual, ΔT)
			Residual_Norm = 0.0::Float64
			 for iZ = 1:Nz
				Residual_Norm += (Residual[iZ] / (ΔT[iT] * discret.ΔZ[iZ])) ^ 2.0
			end # for: iZ=Nz
		return  √(Residual_Norm / Float64(Nz))			
		end  # function: CONVERGENCECRITERIA
	#-----------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SOLVING_TRIAGONAL_MATRIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, hydro, iT::Int64, Nz::Int64, optionHypix, paramHypix, Residual, ΔLnΨmax, θ, Ψ, Ψ_Min::Vector{Float64}, Ψ_Max::Float64)

			Matrix_Trid = Tridiagonal(∂R∂Ψ△[2:Nz], ∂R∂Ψ[1:Nz], ∂R∂Ψ▽[1:Nz-1])

			Residual = reshape(Residual, Nz, 1) # Transforming from row to column

			NewtonStep = Matrix_Trid \ -Residual
			for iZ=1:Nz
				# Iteration k-1
					Ψ₀ = Ψ[iT,iZ]
					θ₀ = θ[iT,iZ]
				
				# Updating Ψ
				if !isnan(NewtonStep[iZ])
					# Newton step
						Ψ[iT,iZ] += NewtonStep[iZ]
				
					# Correction of θ entering a dry soil 
						Ψ = ZHA_WETING_DRYSOIL(hydro, iT, iZ, optionHypix, θ, θ₀, Ψ, Ψ₀)

					# Assuring that the limits of Ψ are physical
						Ψ[iT,iZ] = min(max(Ψ[iT,iZ], Ψ_Min[iZ]), Ψ_Max)

					# Smootening the steps
						Ω = DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT, iZ, optionHypix, ΔLnΨmax, θ₀, Ψ, Ψ₀)
						Ψ[iT,iZ] = Ω * Ψ[iT,iZ] + (1.0 - Ω) * Ψ₀ 
			
				# No comvergence
				else
					# @warn error("===== Difficulties in inverting Tridiagonal =====")
					Ψ[iT,iZ] = Ψ₀ + eps(100.0)
					println(" ================   STRUGGLING ====================")
				end
			end # for iZ=1:Nz	
		return Ψ
		end  # function: SOLVING_TRIAGONAL_MATRIX
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NEWTO_NRAPHSON_STEP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DYNAMIC_NEWTON_RAPHSON_STEP(hydro, iT::Int64, iZ::Int64, optionHypix, ΔLnΨmax, θ₀, Ψ, Ψ₀)

			Ωmax = 1.0 # 1.0
			Ωmin = 0.2 # 0.2

			# Changing sign slow down
			if  sign(Ψ[iT,iZ]) * sign(Ψ₀) == -1
				return Ωmin
			else
				θ₁ = Ψ_2_θDual(optionHypix, Ψ[iT,iZ], iZ, hydro)
				Δθ = abs(θ₁ - θ₀)
				Δθₘₐₓ = timeStep.ΔθMAX(hydro, iT, iZ, optionHypix, ΔLnΨmax, Ψ) 
				
				return Ωmax - (Ωmax - Ωmin)  * min(Δθ / Δθₘₐₓ, 1.0) ^ 2.0
			end
		end  # function: NEWTO_NRAPHSON_STEP
	# ------------------------------------------------------------------
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OVERSHOTTING_WET_DRY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	"""Zha, Y., Yang, J., Yin, L., Zhang, Y., Zeng, W., Shi, L., 2017. A modified Picard iteration scheme for overcoming numerical difficulties of simulating infiltration into dry soil. Journal of Hydrology 551, 56–69. https://doi.org/10.1016/j.jhydrol.2017.05.053 """
			function ZHA_WETING_DRYSOIL(hydro, iT, iZ, optionHypix, θ, θ₀, Ψ, Ψ₀)
				# Ψwet = max( 3.5391 * hydro.σ[iZ]^3 - 20.676 * hydro.σ[iZ]^2 + 24.835 * hydro.σ[iZ] + 15.976, 0.0 )

				Ψwet = max(-2.3116 * hydro.σ[iZ] ^ 2.0 - 2.9372 * hydro.σ[iZ] + 27.83, 0.0)

				Ψdry = exp(1.6216 * log(hydro.σ[iZ]) + 8.7268)

				# Determine if there is any oscilation at the wet or dry end of the θ(Ψ) curve
				# if  (Ψ₀ ≤ Ψwet && Ψ[iT,iZ] ≥ Ψdry) 
				# 	println("Ψ₀ ≤ Ψwet && Ψ[iT,iZ] ≥ Ψdry")
				# end
				if (Ψ[iT,iZ] ≤ Ψwet && Ψ₀ ≥ Ψdry)
					θ[iT,iZ] = θ₀ + (Ψ[iT,iZ] - Ψ₀) * ∂θ∂Ψ(optionHypix, Ψ₀, iZ, hydro)

					θ[iT,iZ] = max(min(θ[iT,iZ], hydro.θs[iZ]), hydro.θr[iZ])

					Ψ[iT,iZ] = θ_2_ΨDual(optionHypix, θ[iT,iZ] , iZ, hydro)
				end  # Ψ[iT,iZ] ≤ Ψwet && Ψ₀ ≥ Ψdry
			return Ψ
			end  # function:ZHA_WETING_DRYSOIL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RERUN_HYPIX
	# 		WITH UPDATED Ψ
	#     Rerun if updated ΔT is smaller compared to previously Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RERUN_HYPIX(discret, Flag_NoConverge::Bool, Hpond::Vector{Float64}, hydro, iCount_ReRun::Int64, iNonConverge::Int64, iT::Int64, Nz::Int64, optionHypix, paramHypix, Pkₐᵥₑᵣ, Q::Matrix{Float64}, ΔLnΨmax::Vector{Float64}, ΔPrThroughfall::Vector{Float64},  ΔPr_Soil, ΔSink::Matrix{Float64}, ΔT::Vector{Float64}, θ::Matrix{Float64}, Ψ::Matrix{Float64})

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function COMPUTE_ΔT(discret, Flag_NoConverge::Bool, Hpond::Vector{Float64}, hydro, iT::Int64, Nz::Int64, optionHypix, paramHypix, Pkₐᵥₑᵣ::Vector{Float64}, Q::Matrix{Float64}, ΔLnΨmax::Vector{Float64}, ΔPr_Soil::Float64, ΔPrThroughfall::Vector{Float64}, ΔSink::Matrix{Float64}, ΔT::Vector{Float64}, θ::Matrix{Float64}, Ψ::Matrix{Float64})

					Q[iT,1] = flux.Q!(optionHypix, discret, hydro, 1, iT, Nz, paramHypix, Pkₐᵥₑᵣ, Hpond, ΔPr_Soil, ΔSink, ΔT, θ, Ψ[iT,1], Ψ[iT,1])
					for iZ=1:Nz
						Q[iT,iZ+1] = flux.Q!(optionHypix, discret, hydro, iZ+1, iT, Nz, paramHypix, Pkₐᵥₑᵣ, Hpond, ΔPr_Soil, ΔSink, ΔT, θ, Ψ[iT, min(iZ+1, Nz)], Ψ[iT,iZ])
					end

					ΔT_New, ~ = timeStep.ADAPTIVE_TIMESTEP(discret, hydro, iT, Nz, optionHypix, paramHypix, Q, ΔLnΨmax, ΔSink, ΔT, θ, Ψ; Flag_NoConverge=Flag_NoConverge)
				return ΔT_New
				end  # function: COMPUTE_ΔT  
			# ------------------------------------------------------------------

			if iCount_ReRun ≤ 3 
				if Flag_NoConverge
				   # We cannot decrease time step further
					if  ΔT[iT] ≤ paramHypix.ΔT_Min + 4.0 
						Flag_ReRun    = false
						iCount_ReRun  = 1
					
					elseif iCount_ReRun == 1
						# Re compute ΔT by taking the smallest ΔT value of the whole prifile
						ΔT[iT] = COMPUTE_ΔT(discret, Flag_NoConverge, Hpond, hydro, iT, Nz, optionHypix, paramHypix, Pkₐᵥₑᵣ, Q, ΔLnΨmax, ΔPr_Soil, ΔPrThroughfall, ΔSink, ΔT, θ, Ψ)
						
                  Flag_ReRun     = true
                  iCount_ReRun   += 1

					# We can decrease the ΔT
					else		
						# ΔT[iT] = paramHypix.ΔT_Min + 0.45 * max(ΔT[iT] - paramHypix.ΔT_Min, 0.0)
						ΔT[iT] = max(0.4 * ΔT[iT], paramHypix.ΔT_Min)
						
						Flag_ReRun     = true
						iCount_ReRun   += 1
					end	
					
				else # Flag_NoConverge
					Δθerror = 0.0
					for  iZ=1:Nz
						ΔθMax = timeStep.ΔθMAX(hydro, iT, iZ, optionHypix, ΔLnΨmax, Ψ)
						Δθerror  += (abs(θ[iT,iZ]- θ[iT-1,iZ]) / ΔθMax) ^ 2.0
					end
					Δθerror  = √(Δθerror / Float64(Nz))

					if Δθerror > 1.0 
						ΔTₒ = COMPUTE_ΔT(discret, Flag_NoConverge, Hpond, hydro, iT, Nz, optionHypix, paramHypix, Pkₐᵥₑᵣ, Q, ΔLnΨmax, ΔPr_Soil, ΔPrThroughfall, ΔSink, ΔT, θ, Ψ)
		
						if ΔTₒ < paramHypix.ΔT_MaxChange * ΔT[iT] 
							Flag_ReRun     = true
							ΔT[iT]         = ΔTₒ
							iCount_ReRun  += 1	
						else # <>=<>=<>=<>=<>
							Flag_ReRun   = false
							iCount_ReRun = 1
						end
					else # Δθerror > 1.0 
						Flag_ReRun   = false
						iCount_ReRun = 1
					end # Δθerror > 1.0 
				end # Flag_NoConverge
			else # if iCount_ReRun ≤ 3
				Flag_ReRun = false
				iCount_ReRun = 1
				if Flag_NoConverge
					iNonConverge += 1
					# println(iNonConverge)
				end
			end  # if iCount_ReRun ≤ 2

		return Flag_ReRun, iCount_ReRun, iNonConverge, ΔT
		end  # function: RERUN_HYPIX
	# ------------------------------------------------------------------

end # module: richard
#......................................................................