# =============================================================
#		MODULE: Δtchange
# =============================================================
module Δtchange
	import ..interpolate, ..tool, ..cst
	import Dates: Second

	function CHANGE_OUTPUT_ΔT(∑Pet, ∑PrThroughfall, ∑T, ∑T_Climate, ∑WaterBalance_η, ∑ΔSink, clim, dateHypix, Hpond, iScenario, Nit::Int64, Nz::Int64, optionHypix, paramHypix, pathInputHypix, Q, ΔEvaporation, ΔRunoff, ΔT, θ, Ψ, ∑ΔQ_Obs, ∑T_Qobs; obsθ=0)

		# PREPROCESSING ∑Evaporation, ∑ΔQ
         ∑Evaporation   = fill(0.0::Float64, Nit)
         ∑Pr_Gross      = fill(0.0::Float64, clim.N_Climate)
         ∑Runoff        = fill(0.0::Float64, Nit)
         ∑ΔQ            = fill(0.0::Float64, Nit, Nz+1)

			∑Evaporation[1]    = ΔEvaporation[1]
         ∑ΔQ[1, 1:Nz+1]  .= 0.0::Float64
			for iT=2:Nit
				∑Evaporation[iT] = ∑Evaporation[iT-1] + ΔEvaporation[iT]
				for iZ=1:Nz+1
					∑ΔQ[iT,iZ] = ∑ΔQ[iT-1, iZ] + ΔT[iT] * Q[iT,iZ]
				end
			end

		# ∑Pr_Gross
			∑Pr_Gross[1] = clim.Pr[1]
			for iT= 2:clim.N_Climate
				∑Pr_Gross[iT] = ∑Pr_Gross[iT-1] + clim.Pr[iT]
			end

		# ∑Runoff
			∑Runoff[1] = ΔRunoff[1]
			for iT=2:Nit
				∑Runoff[iT] = ∑Runoff[iT-1] + ΔRunoff[iT]
			end

		# PREPARING DATA FOR PLOTS
			∑T_Reduced = collect(range(0.0::Float64, step=paramHypix.ΔT_Output, stop=dateHypix.Δ∑T_Sim)) 
			@. ∑T_Reduced = ∑T_Reduced + dateHypix.Δ∑T_StartSim

			Nit_Reduced = length(∑T_Reduced)

		# PREPARING DATES WITH INTERVAL:
			# COMPUTING CUMULATIVE TIME
				Date_Reduced = collect(range(dateHypix.Date_SimStart, step=Second(paramHypix.ΔT_Output), stop=clim.Date[end]))

		# INTERPOLATING DATA
			θ_Reduced = fill(0.0::Float64, Nit_Reduced, Nz)	
				θ_Reduced = interpolate.INTERPOLATE_2D_LOOP(∑T, ∑T_Reduced, Nit, Nz, θ_Reduced, θ)

			if optionHypix.θobs
				θobs_Reduced = fill(0.0::Float64, Nit_Reduced, obsθ.Ndepth)
				θobs_Reduced = interpolate.INTERPOLATE_2D_LOOP(obsθ.∑T[1:obsθ.Nit], ∑T_Reduced, obsθ.Nit, obsθ.Ndepth, θobs_Reduced, obsθ.θobs[1:obsθ.Nit,1:obsθ.Ndepth])
			else
				θobs_Reduced = [0.0]
			end

			Ψ_Reduced = fill(0.0::Float64, Nit_Reduced, Nz)
				Ψ_Reduced =  interpolate.INTERPOLATE_2D_LOOP(∑T, ∑T_Reduced, Nit, Nz, Ψ_Reduced, Ψ) 

			∑ΔQ_Reduced = fill(0.0::Float64, Nit_Reduced, Nz+1)
				∑ΔQ_Reduced = interpolate.INTERPOLATE_2D_LOOP(∑T, ∑T_Reduced, Nit, Nz+1, ∑ΔQ_Reduced, ∑ΔQ)

		# .<>.<>.<>
			ΔT_Reduced = fill(0.0::Float64, Nit_Reduced)
				ΔT_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ΔT_Reduced, ΔT)

			∑Evaporation_Reduced = fill(0.0::Float64, Nit_Reduced)
				∑Evaporation_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ∑Evaporation_Reduced, ∑Evaporation)

			∑Runoff_Reduced = fill(0.0::Float64, Nit_Reduced)
				∑Runoff_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ∑Runoff_Reduced, ∑Runoff)

			∑Pet_Reduced = fill(0.0::Float64, Nit_Reduced)
				∑Pet_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ∑Pet_Reduced, ∑Pet)

			∑PrThroufall_Reduced = fill(0.0::Float64, Nit_Reduced)
				∑PrThroufall_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ∑PrThroufall_Reduced, ∑PrThroughfall)

			∑WaterBalanceη_Reduced = fill(0.0::Float64, Nit_Reduced)
				∑WaterBalanceη_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ∑WaterBalanceη_Reduced, ∑WaterBalance_η)

			∑∑Sink = fill(0.0::Float64, Nit_Reduced)
				∑∑Sink = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, Nit, ∑∑Sink, ∑ΔSink)

			∑PrGross_Reduced = fill(0.0::Float64, Nit_Reduced)	
				∑PrGross_Reduced = interpolate.INTERPOLATE_1D_LOOP(∑T_Climate[1:clim.N_Climate], ∑T_Reduced, Nit_Reduced, clim.N_Climate, ∑PrGross_Reduced, ∑Pr_Gross)

			ΔPond_Reduced = fill(0.0, Nit_Reduced)
				ΔPond_Reduced = interpolate.INTERPOLATE_1D_MAX(∑T, ∑T_Reduced, Nit_Reduced, Nit, ΔPond_Reduced, Hpond)

			if !(isempty(pathInputHypix.Drainage[iScenario]))
				∑ΔQ_Obs_Reduced  = fill(0.0, Nit_Reduced)
				∑ΔQ_Obs_Reduced₀ = fill(0.0, Nit_Reduced)
				∑ΔQ_Obs_Reduced₀ = interpolate.INTERPOLATE_1D_LOOP(∑T_Qobs, ∑T_Reduced, Nit_Reduced, length(∑ΔQ_Obs), ∑ΔQ_Obs_Reduced₀, ∑ΔQ_Obs)

				∑ΔQ_Obs_Reduced[1] = 0.0
				for iT=2:Nit_Reduced
					∑ΔQ_Obs_Reduced[iT] = ∑ΔQ_Obs_Reduced₀[iT] - ∑ΔQ_Obs_Reduced₀[1]
				end
			end

			# From ∑ to Δ
            ΔEvaporation_Reduced   = fill(0.0::Float64, Nit_Reduced)
            ΔQ_Reduced             = fill(0.0::Float64, Nit_Reduced, Nz+1)
            ΔPet_Reduced           = fill(0.0::Float64, Nit_Reduced)
            ΔPrThroughfall_Reduced = fill(0.0::Float64, Nit_Reduced)
            ΔPrGross_Reduced       = fill(0.0::Float64, Nit_Reduced)
            ΔSink_Reduced          = fill(0.0::Float64, Nit_Reduced)
            ΔRunoff_Reduced        = fill(0.0::Float64, Nit_Reduced)

				if !(isempty(pathInputHypix.Drainage[iScenario]))
            	ΔQ_Obs_Reduced         = fill(0.0::Float64, Nit_Reduced)
					ΔQ_Obs_Reduced[1]         = 0.0
				else
					ΔQ_Obs_Reduced = []
				end
				
			# Root Water Uptake daily 
				# Initial condition
            ΔEvaporation_Reduced[1]   = 0.0
            ΔQ_Reduced[1,1:Nz+1]     .= 0.0
            ΔPet_Reduced[1]           = 0.0
            ΔPrThroughfall_Reduced[1] = 0.0
            ΔPrGross_Reduced[1]       = 0.0
            ΔSink_Reduced[1]          = 0.0
            ΔRunoff_Reduced[1]        = 0.0

			# ∑T_Reduced[1] = ∑T_Reduced[1]
			set_zero_subnormals(true)
			@fastmath for iT=2:Nit_Reduced
				for iZ=1:Nz+1
					ΔQ_Reduced[iT,iZ] = ∑ΔQ_Reduced[iT,iZ] - ∑ΔQ_Reduced[iT-1,iZ]
				end

            ΔSink_Reduced[iT]          = ∑∑Sink[iT] - ∑∑Sink[iT-1]

            ΔPet_Reduced[iT]           = ∑Pet_Reduced[iT] - ∑Pet_Reduced[iT-1]

            ΔPrThroughfall_Reduced[iT] = ∑PrThroufall_Reduced[iT] - ∑PrThroufall_Reduced[iT-1]

            ΔPrGross_Reduced[iT]       = max(∑PrGross_Reduced[iT] - ∑PrGross_Reduced[iT-1], 0.0::Float64)
				
            ΔEvaporation_Reduced[iT]   = max(∑Evaporation_Reduced[iT] -  ∑Evaporation_Reduced[iT-1], 0.0::Float64)
				
            ΔRunoff_Reduced[iT]        = max(∑Runoff_Reduced[iT] - ∑Runoff_Reduced[iT-1], 0.0::Float64)

				if !(isempty(pathInputHypix.Drainage[iScenario]))
            	ΔQ_Obs_Reduced[iT]         = max(∑ΔQ_Obs_Reduced[iT] - ∑ΔQ_Obs_Reduced[iT-1], 0.0::Float64)
				else
					∑ΔQ_Obs_Reduced = Float64[]
				end
			end  # for iT=1:Nit
			set_zero_subnormals(false)
	
	return ∑Runoff_Reduced, ∑T_Reduced, ∑WaterBalanceη_Reduced, ∑ΔQ_Obs_Reduced, Date_Reduced, Nit_Reduced, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPrGross_Reduced, ΔPrThroughfall_Reduced, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔSink_Reduced, ΔT_Reduced, θ_Reduced, θobs_Reduced, Ψ_Reduced
	end # function: CHANGE_OUTPUT_ΔT
	
end  # module: Δtchange
# ............................................................