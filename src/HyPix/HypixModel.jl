# =============================================================
#		module: hypix
# =============================================================
module hypixModel

	import ..evaporation, ..interception, ..interpolate, ..pet, ..richard, ..rootWaterUptake, ..sorptivity, ..timeStep, ..ΨminΨmax
	import ..wrc: θ_2_ΨDual, Ψ_2_θDual

	export HYPIX_MODEL

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIX_MODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYPIX_MODEL(∂K∂Ψ::Vector{Float64}, ∂R∂Ψ::Vector{Float64}, ∂R∂Ψ△::Vector{Float64}, ∂R∂Ψ▽::Vector{Float64}, ∑Pet_Climate::Vector{Float64}, ∑Pet::Vector{Float64}, ∑PrThroughfall_Climate::Vector{Float64}, ∑PrThroughfall::Vector{Float64}, ∑T_Climate::Vector{Float64}, ∑T::Vector{Float64}, clim, CropCoeficientᵀ_η::Vector{Float64}, CropCoeficientᵀ::Vector{Float64}, discret, Flag_θΨini::Symbol, Hpond::Vector{Float64}, hydro, iScenario::Int64, K_Aver_Vect::Vector{Float64}, K_Aver₀_Vect::Vector{Float64}, N_∑T_Climate::Int64, N_iRoot::Int64, Nz::Int64, optionHypix, paramHypix, pathInputHypix, Pkₐᵥₑᵣ::Vector{Float64}, Q::Matrix{Float64}, Residual::Vector{Float64}, veg, Z::Vector{Float64}, ΔEvaporation::Vector{Float64}, ΔLnΨmax::Vector{Float64}, ΔPet::Vector{Float64}, ΔPrThroughfall::Vector{Float64}, ΔRunoff::Vector{Float64}, ΔSink::Matrix{Float64}, ΔRootDensity, ΔT::Vector{Float64}, θ::Matrix{Float64}, θini_or_Ψini::Vector{Float64}, Ψ_Min::Vector{Float64}, Ψ::Matrix{Float64}, Ψbest::Vector{Float64})

		# VEGETATION PARAMETERS WHICH VARY WITH TIME

			for iT = 1:clim.N_Climate
				# if optionHypix.LookupTable_Lai
				# 	Laiᵀ[iT] = (veg.Lai_Max - veg.Lai_Min) * Laiᵀ_η[iT] + veg.Lai_Min
				# else
				# 	Laiᵀ[iT] = veg.Lai
				# end
				if optionHypix.LookUpTable_CropCoeficient
					CropCoeficientᵀ[iT]  = (veg.CropCoeficient_Max - veg.CropCoeficient_Min) * CropCoeficientᵀ_η[iT] + veg.CropCoeficient_Min
				else
					CropCoeficientᵀ[iT]  = veg.CropCoeficient
				end
			end # for

		# RAINFALL INTERCEPTION
			if optionHypix.RainfallInterception
				∑Pet_Climate, ∑PrThroughfall_Climate, clim = interception.RAINFALL_INTERCEPTION_START(∑Pet_Climate, ∑PrThroughfall_Climate, clim, optionHypix, veg)
			end
		
		# ROOTS
		if optionHypix.RootWaterUptake
			# If root density function is computed by model, or else the root density is read from file
			if isempty(pathInputHypix.∑RootDensityFunc[iScenario])
				N_iRoot = rootWaterUptake.rootDistribution.N_IROOT(Nz, veg, Z) # Last cell of rootzone
				ΔRootDensity = rootWaterUptake.rootDistribution.ROOT_DENSITY(discret, N_iRoot, veg, Z)
			end
		else
			ΔRootDensity = 0.0::Float64
			N_iRoot = 1::Int64
		end # optionHypix.RootWaterUptake

		# if optionHypix.Evaporation 
		# 	N_iEvapo = evaporation.N_IEVAPO(Nz, veg, Z) # Smap_Depth where evaporation can occure
		# end # optionHypix.Evaporation

		# MINIMUM OR MAXIMUM Ψ VALUES THIS IS SUCH THAT ∂Θ∂Ψ ≠ 0 WHICH INFLUENCES THE NEWTON-RAPHSON METHOD TO BE REMOVED
			Ψ_Max = 0.0
			for iZ=1:Nz
            Ψ_Max₀,~ = ΨminΨmax.ΨMINΨMAX(hydro.θs[iZ], hydro.θsMacMat[iZ],  hydro.σ[iZ],  hydro.σMac[iZ], hydro.Ψm[iZ], hydro.ΨmMac[iZ])
				Ψ_Max = max(Ψ_Max, Ψ_Max₀)
			end  # for iZ=1:Nz

		# ADAPTIVETIMESTEP
			ΔLnΨmax = timeStep.ΔΨMAX!(hydro, Nz, optionHypix, paramHypix, ΔLnΨmax)

		# FIRST TIME STEP
         Flag_NoConverge         = false::Bool
         Flag_ReRun              = false::Bool
         IterCount               = 0::Int64
         iTer                    = 10
         iNonConverge            = 0::Int64
         iT                      = 1::Int64
         iT_Pet                  = 2::Int64
         iT_Pr                   = 2::Int64
         ΔEvaporation[1]         = 0.0::Float64
         Hpond[1]                = 0.0::Float64
         ΔPet[1]                 = 0.0::Float64
         ΔPrThroughfall[1]       = 0.0::Float64
         ΔSink[1,1:Nz]           .= 0.0::Float64
         ΔT[1]                   = 0.0::Float64
         ∑Pet[1]                 = 0.0::Float64
         ∑PrThroughfall[1]       = 0.0::Float64
         ∑T[1]                   = 0.0::Float64
         iCount_ReRun            = 1::Int64
			
		# Boundary conditions
			if Flag_θΨini == :θini
				for iZ = 1:Nz
               θ[1,iZ] = max( min(hydro.θs[iZ] * 0.9, θini_or_Ψini[iZ]), min(hydro.θr[iZ] * 1.1, hydro.θs[iZ] * 0.9) ) # Just in case
               Ψ[1,iZ] = θ_2_ΨDual(optionHypix, θ[1,iZ], iZ, hydro)
				end

			elseif Flag_θΨini == :Ψini
				for iZ = 1:Nz
               Ψ[1,iZ] = θini_or_Ψini[iZ]
               θ[1,iZ] = Ψ_2_θDual(optionHypix, Ψ[1,iZ], iZ, hydro)
				end
			end

			if optionHypix.TopBoundary⍰ == "Ψ"
				Ψ[1,1] = paramHypix.Ψ_Top
				θ[1,1]  = Ψ_2_θDual(optionHypix, Ψ[1,1], 1, hydro)
			end

			if optionHypix.BottomBoundary⍰ == "Ψ"
				Ψ[1,Nz] = paramHypix.Ψ_Botom
				θ[1,Nz]  = Ψ_2_θDual(optionHypix, Ψ[1,Nz], Nz, hydro)	
			end

			for iZ = 1:Nz
            Ψbest[iZ] = Ψ[1,iZ]
            Q[1,Nz]  = 0.0::Float64
			end
			Q[1,Nz+1] = 0.0::Float64

		# =+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+
		while true # this controles the time loop
			# INCREASING OR DECREASING THE TIME STEP
				∑T, FlagContinueLoop, iT, ΔT, Δθ_Max = timeStep.TIMESTEP(∑T, discret, Flag_ReRun, hydro, iT, iTer, Float64(N_∑T_Climate), Nz, optionHypix, paramHypix, Q, ΔLnΨmax, ΔSink, ΔT, θ, Ψ)

				if FlagContinueLoop == false
					iT = iT - 1
					break # End of simulation
				end

			# DERIVING FORCING DATA ΔPrThroughfall & ΔPet:
				∑PrThroughfall[iT], ΔPrThroughfall[iT], iT_Pr = interpolate.∑_2_Δ(∑PrThroughfall[iT-1], ∑PrThroughfall_Climate, ∑T, ∑T_Climate, iT_Pr, clim.N_Climate, Flag_ReRun, iT)

			# POTENTIAL EVAPOTRANSPIRATION
				if optionHypix.RootWaterUptake || optionHypix.Evaporation
					∑Pet[iT], ΔPet[iT], iT_Pet = interpolate.∑_2_Δ(∑Pet[iT-1], ∑Pet_Climate, ∑T, ∑T_Climate, iT_Pet, clim.N_Climate, Flag_ReRun, iT)
				end # optionHypix.RootWaterUptake || optionHypix.Evaporation

				if optionHypix.Evaporation						
					ΔPet_Evap, ΔPet_Transp = pet.BEER_LAMBERT_LAW(iT, iT_Pr, ΔPet, veg, clim)
					
				else
					ΔPet_Transp = ΔPet[iT]
					ΔPet_Evap = 0.0::Float64
				end
				
			# ROOT WATER UPTAKE MODEL
				if optionHypix.RootWaterUptake
					ΔSink = rootWaterUptake.ROOT_WATER_UPTAKE(CropCoeficientᵀ[iT_Pr-1], iT, N_iRoot, optionHypix, veg, ΔPet_Transp, ΔRootDensity, ΔSink, Ψ)					
				end # optionHypix.RootWaterUptake

			# EVAPORATION FROM THE SURFACE WITH HIGHEST Se
				if optionHypix.Evaporation
					ΔEvaporation = evaporation.EVAPORATION!(hydro, iT, ΔEvaporation, ΔPet_Evap, θ)
					
					ΔSink[iT,1] += ΔEvaporation[iT]
				end # optionHypix.Evaporation

			# Checking that not too much water is removed from the layer
				if optionHypix.RootWaterUptake || optionHypix.Evaporation
					for iZ=1:N_iRoot
						ΔSink[iT,iZ] = min(ΔSink[iT,iZ], discret.ΔZ[iZ] * (θ[iT-1,iZ] - hydro.θr[iZ]))
					end
				end # if: optionHypix

			# SORPTIVITY TO COMPUTE INFILTRATION RATE
				Sorptivity = sorptivity.SORPTIVITY(θ[iT-1, 1], 1, hydro, optionHypix, optionHypix; Rtol = 10^-3.0, SorptivityModelScaled=false)

			# SOLVING THE EXPLICIT RICHARDS
			Flag_NoConverge, Flag_ReRun, Hpond, iCount_ReRun, iNonConverge, iTer, IterCount, Q, ΔRunoff, ΔT, θ, Ψ = richard.RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, Hpond, hydro, iCount_ReRun, iNonConverge, iT, IterCount, K_Aver_Vect, K_Aver₀_Vect, Nz, optionHypix, paramHypix, Pkₐᵥₑᵣ, Q, Residual, Sorptivity, ΔLnΨmax, ΔPrThroughfall, ΔRunoff, ΔSink, ΔT, θ, Ψ_Max, Ψ_Min, Ψ, Ψbest)
				
			# SPECIAL BOUNDARY CONDITIONS
				if optionHypix.TopBoundary⍰ == "Ψ"
					ΔPrThroughfall[iT] = ΔT[iT] * Q[iT, 1]
					∑PrThroughfall[iT] = ∑PrThroughfall[iT-1] + ΔPrThroughfall[iT]
				end
		end # while loop
		# =+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+=+=+=+==+=+=+=+=+=+	

		Nit = iT # Maximum time steps

	return ∑Pet, ∑PrThroughfall, ∑T, ∑T_Climate, clim, discret, Hpond, iNonConverge, IterCount, N_iRoot, Nit, Nz, Q, veg, ΔEvaporation, ΔRootDensity, ΔRunoff, ΔT, θ, Ψ
	end  # function: HYPIX_MODEL
	
end  # module hypix
# ............................................................