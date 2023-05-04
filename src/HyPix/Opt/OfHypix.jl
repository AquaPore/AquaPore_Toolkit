# =============================================================
#		MODULE: ofHypix
# =============================================================
module ofHypix

		import ..tool, ..interpolate, ..stats
		import Statistics
		export WOF_θ, RMSE_θ, CCC_θ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : WOF_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function WOF_θ(∑T, Nit::Int, Nz::Int, obsθ, paramHypix, Hpond, θ, θSim; 🎏_WofDepth=false)

				θSim = interpolate.INTERPOLATE_2D_LOOP(∑T, obsθ.∑T[1:obsθ.Nit], Nit, Nz, θSim, θ)
				
				Wof = 0.0
				iCount = 0

				for iZ=1:obsθ.Ndepth

					if 🎏_WofDepth
						Wdepth = 2.0 * (Float64(obsθ.Ndepth) + 1.0 - Float64(iZ) ) / (Float64(obsθ.Ndepth) * (Float64(obsθ.Ndepth) + 1.0))
					else
						Wdepth = 1.0
					end

					for iT=1:obsθ.Nit		
						if θSim[iT,iZ] > 0.0 && obsθ.θobs[iT,iZ] > 0.0# avoiding no data

							θerror = abs(obsθ.θobs[iT,iZ] - θSim[iT, obsθ.ithetaObs[iZ]])
		
							Wof += Wdepth * θerror ^ 2.0
							# end # if θobs_Min ≤ θSim[iT, obsθ.ithetaObs[iZ]] ≤ θobs_Max

							iCount += 1

						end # if: obsθ.θobs[iT,iZ] > 0.0
					end # for iT
				end # for iZ

				# Penalty if we have too much ponding
				Wof_Pond = max(Hpond[Nz] - paramHypix.opt.ΔHpondMax, 0.0) / 1000.0

			return Wof = √(Wof / Float64(iCount)) + Wof_Pond
			end # function WOF_θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : RMSE_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function RMSE_θ(∑T, obsθ, Nit::Int, Nz::Int, θ, θSim)
				θSim = interpolate.INTERPOLATE_2D_LOOP(∑T, obsθ.∑T[1:obsθ.Nit], Nit, Nz, θSim, θ)
				
				Rmse = 0.0
				iCount = 0
				for iZ=1:obsθ.Ndepth
					for iT=1:obsθ.Nit 	
						if θSim[iT,iZ] > 0.0 && obsθ.θobs[iT,iZ] > 0.0 # avoiding no data
							Rmse += (obsθ.θobs[iT,iZ] - θSim[iT, obsθ.ithetaObs[iZ]]) ^ 2.0
							iCount += 1
						end # if: obsθ.θobs[iT,iZ] > 0.0
					end # for it

				end # for iZ

			return √(Rmse / (Float64(iCount)))		
			end # function WOF_θ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : OF_Q
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			""" Compares Hypix and simulated cumulative recharge at time step ΔTq"""
			function OF_Q(∑T::Vector{Float64}, ∑T_Qobs::Vector{Float64}, ∑ΔQ_Obs::Vector{Float64}, Nit::Int64, Nz::Int64, obsθ, Q::Array{Float64}, ΔT::Vector{Float64}; ΔTq=60*60*24*30)

				∑T_StartCalibration = obsθ.∑T[1]

				# NEW TIME STEP WHERE WE ARE COMPUTTING CUMULATIVE Q
					# Number of time steps ΔTq
						Nit_Reduced = Int64((∑T[Nit] - ∑T_StartCalibration) ÷ ΔTq)

					# Time where we comute Q
						∑T_Reduced = collect(∑T_StartCalibration: ΔTq : Nit_Reduced * ΔTq + ∑T_StartCalibration)

						Nit_Reduced = length(∑T_Reduced) # Updated since there could be difference in 1
						
				# CUMULATIVE HYPIX ΔQ
					∑ΔQbottom = fill(0.0::Float64, Nit)
					∑ΔQbottom[1] = ΔT[1] * Q[1,Nz+1]
						for iT=2:Nit
							∑ΔQbottom[iT] = ∑ΔQbottom[iT-1] + ΔT[iT] * Q[iT,Nz+1]
						end  # for iT=1:Nit

				# INTERPOLATION 
					# Observed data
						∑ΔQ_Obs_Interp = fill(0.0::Float64, Nit_Reduced)
						∑ΔQ_Obs_Interp = interpolate.INTERPOLATE_1D_LOOP(∑T_Qobs, ∑T_Reduced, Nit_Reduced, length(∑T_Qobs), ∑ΔQ_Obs_Interp, ∑ΔQ_Obs)

					# HyPix data
						∑ΔQ_Hypix_Interp = fill(0.0::Float64, Nit_Reduced)
						∑ΔQ_Hypix_Interp = interpolate.INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, Nit_Reduced, length(∑T), ∑ΔQ_Hypix_Interp, ∑ΔQbottom)

				# STATISTICS
               ΔQ_Obs_Interp      = fill(0.0, Nit_Reduced)
               ΔQ_Hypix_Interp    = fill(0.0, Nit_Reduced)
               ΔQ_Obs_Interp[1]   = ∑ΔQ_Obs_Interp[1]
               ΔQ_Hypix_Interp[1] = ∑ΔQ_Hypix_Interp[1]
					for iT=2:Nit_Reduced
                  ΔQ_Obs_Interp[iT]   = ∑ΔQ_Obs_Interp[iT] - ∑ΔQ_Obs_Interp[iT-1]
                  ΔQ_Hypix_Interp[iT] = ∑ΔQ_Hypix_Interp[iT] -  ∑ΔQ_Hypix_Interp[iT-1]
					end 

					Of_Q = max(stats.NSE_MINIMIZE(∑ΔQ_Hypix_Interp, ∑ΔQ_Obs_Interp; Power=2.0), 0.0)
				
			return Of_Q
			end  # function: OF_Q
		# ------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : CCC_θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function CCC_θ(∑T, Nit::Int, Nz::Int, obsθ, θ, θSim)

				θSim = interpolate.INTERPOLATE_2D_LOOP(∑T, obsθ.∑T[1:obsθ.Nit], Nit, Nz, θSim, θ)

				θobs_Layer = fill(0.0::Float64, obsθ.Nit)
				θSim_Layer = fill(0.0::Float64, obsθ.Nit)

				WofCcc = 0.0

				for iZ=1:obsθ.Ndepth
					iTgood = 1
					for iT=1:obsθ.Nit
						if θSim[iT,iZ] > 0.0 && obsθ.θobs[iT,iZ] > 0.0 # avoiding no data
							θobs_Layer[iTgood] = obsθ.θobs[iT,iZ]
							θSim_Layer[iTgood] = θSim[iT,obsθ.ithetaObs[iZ]]
							iTgood += 1
						end # if: obsθ.θobs[iT,iZ] > 0.0
					end # for iT

					iTgood -= 1

					WofCcc += stats.RMSE_CONCORDANCE_CORELATION_COEFICIENT(θobs_Layer[1:iTgood], θSim_Layer[1:iTgood])
				end # for iZ

			return WofCcc
			end # function WOF_θ

end  # module ofHypix
# ............................................................