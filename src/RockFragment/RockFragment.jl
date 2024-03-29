# =============================================================
#		module: rockFragment
# =============================================================
module rockFragment
	import Polynomials

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  ρᵦ_2_Φ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρᵦ_2_Φ(hydro, NiZ::Int64, option, RockFragment::Vector{Float64}, ρₚ_Fine::Vector{Float64}, ρₚ_Rock::Vector{Float64}, ρᵦ_Soil::Vector{Float64})

			if option.rockFragment.RockInjectedIncluded⍰  == "InjectRock"
				return hydro = rockFragment.injectRock.ρᵦ_2_Φ(hydro, NiZ, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)

			elseif option.rockFragment.RockInjectedIncluded⍰  == "Included"
				return hydro = rockFragment.included.ρᵦ_2_Φ(hydro, NiZ, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)
			end
		end  # function: function ρᵦ_2_Φ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CORECTION_θΨ_WETABLE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CORECTION_θΨ_WETABLE!(N_θΨobs::Vector{Int64}, NiZ::Int64, OrganicMatter::Vector{Float64}, rfWetable, RockClass::Vector{String}, RockFragment::Vector{Float64}, θ_θΨobs::Matrix{Float64}, Ψ_θΨobs::Matrix{Float64})
			for iZ = 1:NiZ	
				iRockClass = rfWetable.RockClass_Dict[RockClass[iZ]]

				WettableFunction = Polynomials.Polynomial(rfWetable.RockClass_Polynomial_Array[iRockClass][1])

				for iθ=1:N_θΨobs[iZ]
					θ_Wettable = WettableFunction(log1p(Ψ_θΨobs[iZ,iθ]))

					# Special case for Pumice soils
					if θ_Wettable < 0.0
						θ_Wettable = θ_θΨobs[iZ,iθ] * 0.75 
					end

					θ_θΨobs[iZ,iθ] = θ_θΨobs[iZ,iθ] + θ_Wettable * RockFragment[iZ]
				end # for iθ=1:N_θΨobs[iZ]
			end #  for iZ = 1:NiZ			
		return  θ_θΨobs
		end  # function: CORECTION_θΨ_WETABLE
	

	# =============================================================
	#		module: injectRock
	# =============================================================
	module injectRock
		import Polynomials
		export CORECTION_Φ!, CORECTION_θΨ!, ρᵦ_2_Φ, CORECTION_Φ!, ϜUNC_RF_CORECTION_KΨ_NASEI, RF_CORECTION_KΨ!

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : CORECTION_θΨ!
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function CORECTION_θΨ!(N_θΨobs::Vector{Int64}, NiZ::Int64, RockFragment::Vector{Float64}, θ_θΨobs::Array{Float64})
				for iZ = 1:NiZ 
					for iθ=1:N_θΨobs[iZ]
						θ_θΨobs[iZ,iθ] = θ_θΨobs[iZ,iθ] * (1.0 - RockFragment[iZ])
					end # for iθ=1:N_θΨobs[iZ]
				end #  for iZ = 1:NiZ	
			return θ_θΨobs
			end  # function: STONECORRECTION_NONEWETABLE
		#..................................................................


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : CORECTION_KΨ!
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function RF_CORECTION_KΨ!(NiZ::Int64, N_KΨobs::Vector{Float64}, RockFragment::Vector{Float64}, K_KΨobs::Vector{Float64}; RF_Start_Increase=0.35, RF_End_Increase=0.8)

				for iZ=1:NiZ, iΨ=1:N_KΨobs[iZ]
					if RockFragment[iZ] < RF_Start_Increase

						K_KΨobs[iZ,iΨ] = K_KΨobs[iZ,iΨ] * ϜUNC_RF_CORECTION_KΨ_NASEI(RockFragment[iZ])
					
					# Ks starts to increase when Ks > RF_Start_Increase linearly, nevertheless due to that there is no data available to test the model, this is an empirical model which has not been validated/tested.
					else
						X = [RF_Start_Increase, RF_End_Increase]
						Y = [ϜUNC_RF_CORECTION_KΨ_NASEI(RF_Start_Increase), 1.0]

						 ϜUNC_INTERPOLATE = Polynomials.fit(X, Y, 1)

						K_KΨobs[iZ,iΨ] = K_KΨobs[iZ,iΨ] * ϜUNC_RF_CORECTION_KΨ_NASEI(RF_Start_Increase) * ϜUNC_INTERPOLATE(min(RockFragment[iZ], RF_End_Increase))
					end
				end #  for	
			return  K_KΨobs
			end  # function: STONECORRECTION_NONEWETABLE
		#..................................................................


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : FUNC_RF_CORECTION_KΨ_NASEI
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			""" Naseri (2022) Rock fragments influence the water retention and hydraulic conductivity of soils """
			function ϜUNC_RF_CORECTION_KΨ_NASEI(RockFragment::Float64; Pshape=1.26, Pcritical=1.26)
				return (1.0 - RockFragment / Pcritical) ^ Pshape	
			end  #  FUNC_RF_CORECTION_KΨ_NASEI()
		# ------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION :  CORECTION_Φ!
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function CORECTION_Φ!(hydro, NiZ, option, RockFragment)		
				if option.run.RockCorection
					for iZ = 1:NiZ 
						hydro.Φ[iZ] = hydro.Φ[iZ] * (1.0 - RockFragment[iZ])
					end
				end
			return  hydro
			end  # function: STONECORRECTION_NONEWETABLE
		#..................................................................


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION :  ρᵦ_2_Φ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ρᵦ_2_Φ(hydro, NiZ::Int64, option, RockFragment::Vector{Float64}, ρₚ_Fine::Vector{Float64}, ρₚ_Rock::Vector{Float64}, ρᵦ_Soil::Vector{Float64})

				for iZ=1:NiZ
					if option.run.RockCorection
						hydro.Φ[iZ] = (1.0 - ρᵦ_Soil[iZ] / ρₚ_Fine[iZ]) * (1.0 - RockFragment[iZ])

					else
						hydro.Φ[iZ] = 1.0 - ρᵦ_Soil[iZ] / ρₚ_Fine[iZ]
					
					end
				end # for
			return hydro
			end  # function: Φ			
		end  # module: injectRock
		# ............................................................


		# =============================================================
		#		module: included
		# =============================================================
		module included
		export ρᵦ_2_Φ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION :  ρᵦ_2_Φ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ρᵦ_2_Φ(hydro, NiZ, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)

				for iZ=1:NiZ
					if option.run.RockCorection					
						hydro.Φ[iZ] = 1.0 - (RockFragment[iZ] * ρᵦ_Soil[iZ] / ρₚ_Rock[iZ]) - ((1.0 - RockFragment[iZ]) * ρᵦ_Soil[iZ] / ρₚ_Fine[iZ])
				
					else
						hydro.Φ[iZ] = 1.0 - ρᵦ_Soil[iZ] / ρₚ_Fine[iZ]
					end
				end # for
			return hydro
			end  # function: Φ
		
	end  # module: included
	# ............................................................
	
end  # module: rockFragment
# ............................................................