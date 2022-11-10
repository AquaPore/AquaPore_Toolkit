	# =============================================================
	#		module non Core : smap
	# =============================================================
	module readSmap
	   import ..tool
   	import Polynomials
		using CSV, Tables, DataFrames
		export SMAP, ROCKFRAGMENT_WETTABLE_STRUCT, IMPERMEABLE_CLASS

			const RockFragment_Max = 0.9

			function SMAP(IdSelect, NiZ, path)
				println("    ~  $(path.inputSmap.Smap) ~")

				Data = CSV.read(path.inputSmap.Smap, DataFrame, header=true)

				DataFrames.sort!(Data, [:Id])
				
            IsTopsoil              = convert(Vector{Float64}, Data."IsTopsoil")
            IsTopsoil              = Int64.(IsTopsoil)
            IsTopsoil              = IsTopsoil[IdSelect]

            Soilname               = convert(Vector{String}, Data."Soilname")
            Soilname               = Soilname[IdSelect]
			
            RockClass              = convert(Vector{String}, Data."RockClass")
            RockClass              = RockClass[IdSelect]
				
            Smap_Depth             = convert(Vector{Float64}, Data."depth_mm")
            Smap_Depth             = Smap_Depth[IdSelect]

            RockFragment           = convert(Vector{Float64}, Data."Stone_Prop")
            RockFragment           = min.(RockFragment_Max, RockFragment)
            RockFragment           = RockFragment[IdSelect]

            Smap_RockDepth         = convert(Vector{Float64}, Data."RockDepth_mm")
            Smap_RockDepth         = Smap_RockDepth[IdSelect]
				
            Smap_MaxRootingDepth   = convert(Vector{Float64}, Data."MaxRootingDepth_mm")
            Smap_MaxRootingDepth   = Smap_MaxRootingDepth[IdSelect]

            Smap_SmapFH            = convert(Vector{String}, Data."SmapFH")
            Smap_SmapFH            = Smap_SmapFH[IdSelect]

            Smap_PermeabilityClass = convert(Vector{String}, Data."PermeabilityClass")
            Smap_PermeabilityClass = Smap_PermeabilityClass[IdSelect]

            Smap_IsImpermeable     = convert(Vector{Bool}, Data."IsImpermeable")
            Smap_IsImpermeable     = Smap_IsImpermeable[IdSelect]

            Smap_ImpermClass       = convert(Vector{String}, Data."ImpermClass")
            Smap_ImpermClass       = Smap_ImpermClass[IdSelect]

				NiZ = length(Smap_ImpermClass)

				# KS_IMPERMEABLE
					KsImpClass_Dict = readSmap.IMPERMEABLE_CLASS(path.inputSmap.LookupTable_Impermeable)

					Ks_Impermeable = fill(-1.0, NiZ)
					for iZ=1:NiZ
						if Smap_IsImpermeable[iZ]
							Ks_Impermeable[iZ] = KsImpClass_Dict[Smap_ImpermClass[iZ]]
						end
					end
			
			return IsTopsoil, Ks_Impermeable, RockClass, RockFragment, Smap_Depth, Smap_MaxRootingDepth, Smap_PermeabilityClass, Smap_RockDepth, Smap_SmapFH, Soilname
			end  # function: SMAP


		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# #		FUNCTION :SMAP
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 	struct SMAP_STRUCT
		# 		Smap_Depth        ::Vector{Float64}
		# 		IsTopsoil    ::Vector{Int64}
		# 		Soilname     ::Vector{String}
		# 		RockFragment ::Vector{Float64}
		# 		RockClass    ::Vector{String}
		# 		Smap_RockDepth    ::Vector{Float64}
		# 		Smap_MaxRootingDepth ::Vector{Float64}
		# 	end

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ROCKFRAGMENT_WETTABLE
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			struct ROCKFRAGMENT_WETTABLE_STRUCT
				# RockClass::Array{String}
				RockClass_Dict::Dict{String, Int64} 
				Ψ_Rf::Array{Float64} 
				θ_Rf::Array{Float64}
				N_Ψ::Array{Int64}
				N_RockClass::Int64
				RockClass_Polynomial_Array::Array{} 
			end
			function ROCKFRAGMENT_WETTABLE(Path)
				println("    ~  $(Path) ~")
				
				# Read data
					Data = CSV.read(Path,  DataFrame, header=true)

					RockClass = convert(Vector{String}, Data. "RockClass")

					N_RockClass = length(RockClass)

					RockClass_Unique = unique(RockClass)
					
					N_RockClass = length(RockClass_Unique)

				# Dictionary
					RockClass_Dict = Dict("a"=>9999)
					for i=1:N_RockClass
						RockClass_Dict[RockClass_Unique[i]] = i
					end

				# Read data
					Ψ₂ = convert(Vector{Float64}, Data."H[mm]")
					θ₂ = convert(Vector{Float64}, Data."Theta[0-1]") 

					N₂ =length(Ψ₂)

					Ψ_Rf = zeros(Int64,(N_RockClass, 100))
					θ_Rf = zeros(Float64,(N_RockClass, 100))
					N_Ψ = zeros(Int64,(N_RockClass))

					iRockClass=1 ; iΨ=1
					for i=1:N₂
						if RockClass[i] == RockClass_Unique[iRockClass]
							Ψ_Rf[iRockClass,iΨ] = Ψ₂[i]
							θ_Rf[iRockClass,iΨ] = θ₂[i]
						else
							N_Ψ[iRockClass]  = iΨ -1
							iRockClass += 1
							iΨ = 1
							Ψ_Rf[iRockClass,iΨ] = Ψ₂[i]
							θ_Rf[iRockClass,iΨ] = θ₂[i]
						end
						iΨ += 1
					end # for i=1:N

					N_Ψ[iRockClass]  = iΨ - 1

				RockClass_Polynomial_Array = []
				for iRockClass=1:N_RockClass
					RockClass_Polynomial = Polynomials.fit(log1p.(Ψ_Rf[iRockClass,1:N_Ψ[iRockClass]]), θ_Rf[iRockClass,1:N_Ψ[iRockClass]])
					X = log1p.(Ψ_Rf[iRockClass,1:N_Ψ[iRockClass]])

					Coeffs = Polynomials.coeffs(RockClass_Polynomial)
				
					RockClass_Polynomial_Array = push!(RockClass_Polynomial_Array, [Coeffs])
				end

			return ROCKFRAGMENT_WETTABLE_STRUCT(RockClass_Dict, Ψ_Rf, θ_Rf, N_Ψ, N_RockClass, RockClass_Polynomial_Array)	
			end  # function: ROCKFRAGMENT_WETTABLE
		#------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Impermeable class
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function IMPERMEABLE_CLASS(Path)
					
				# Read data
					Data = CSV.read(Path,  DataFrame, header=true)

               ImpermClass    = convert(Vector{String}, Data."ImpermClass")
               Ks_ImpermClass = convert(Vector{Float64}, Data."Ks[mm s-1]")
               N_ImpermClass  = length(ImpermClass)

				# Dictionary
					KsImpClass_Dict = Dict("a"=>0.1)
					for i=1:N_ImpermClass
						KsImpClass_Dict[ImpermClass[i]] = Ks_ImpermClass[i]
					end

			return KsImpClass_Dict
			end
		#------------------------------------------------------------------
		
	end  # module: smap
	# ...........................................................