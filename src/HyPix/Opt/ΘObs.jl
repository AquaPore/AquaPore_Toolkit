# =============================================================
#		module: thetaObs
# =============================================================
module thetaObs
	import Dates: value, DateTime
	export thetaObs
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΘOBS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ΘOBS(obsθ, clim, discret, Z::Vector{Float64})

			# CHECKING DATA CONSISTENCY
				if obsθ.Date[1] < clim.Date[2]
					error("\n Hypix error: Starting date of obsθ  $(obsθ.Date[1]) < starting date of climate data $(clim.Date[1]) by 2 iT")
				end # Error checking

				# Checking the celLs
				if obsθ.Z[obsθ.Ndepth] > Z[discret.Nz]
					error("\n Hypix error: depth of measured θ deeper than the max depth of discretisation: obsθ.Z[obsθ.Ndepth]= $(obsθ.Z[obsθ.Ndepth]) > discret.Z[discret.Nz] =$(discret.Z[discret.Nz])") 
				end

			# COMPUTING CUMULATIVE TIME
				for iT=1:obsθ.Nit
					obsθ.∑T[iT] = value(obsθ.Date[iT] - clim.Date[1]) / 1000
				end  # for it=1:obsθ.Nit

				# TRANSFORM THE DEPTH OF MEASURED Θ -> CELL DEPTH
				for iDepth = 1:obsθ.Ndepth
					for iZ = 1:discret.Nz
						if iZ == 1
							if 0.0 ≤ obsθ.Z[iDepth] ≤ Z[1]
								obsθ.ithetaObs[iDepth] = 1
								break  
							end  # if: discret.Z_CellUp
						elseif iZ ≠ 1
							if Z[iZ-1] ≤ obsθ.Z[iDepth] ≤ Z[iZ]
								obsθ.ithetaObs[iDepth] = iZ
								break  
							end  # if: discret.Z_CellUp
						end # if iZ == 1
					end # iZ = 2:discret.Nz						
				end  # for iDepth = 1:obsθ.Ndepth
		
		return obsθ
		end  # function: θOBS
	# .........................................................................................
end # module: thetaObs