# =============================================================
#		module horizonLayer
# =============================================================
module horizonLayer

	import ..hydroStruct
	export HYDROHORIZON_2_HYDRO

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :   HORIZON_2_LAYER
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDROHORIZON_2_HYDRO(hydroHorizon₀, Layer₀::Vector{Float64}, Nz::Int64, optionHypix)

			# Making the numbering of the layers continous without the 0.5 for smootening kayers
			
			Layer_Int = fill(0::Int64, Nz)
			
			iCount = 1
			Layer_Int[1] = 1
			for iZ=2:Nz
				if Layer₀[iZ] > Layer₀[iZ-1]
					iCount += 1
				end
				Layer_Int[iZ] = iCount
			end

			hydro = hydroStruct.HYDROSTRUCT(optionHypix, Nz)

			# Field names of the structure
				FieldName_Array = propertynames(hydroHorizon₀)

				length(hydroHorizon₀.θs)
			
			# looping through every fieldnames of the structure
				for FieldName in FieldName_Array
					Vector = fill(0.0::Float64, Nz)
					for iZ = 1:Nz
						Vector[iZ] = getfield(hydroHorizon₀, FieldName)[Layer_Int[iZ]]
					end
					setfield!(hydro, Symbol(FieldName), Vector)
				end
		return hydro
		end # function HORIZON_2_LAYER
	#---------------------------------------------------------------------

end  # mmodule horizonLayer
# ............................................................