	# =============================================================
	#		module non Core : smap
	# =============================================================
	module readSmap

	   import ..tool
   	import Polynomials
		using CSV, Tables, DataFrames
		export SMAP, ROCKFRAGMENT_WETTABLE_STRUCT, IMPERMEABLE_CLASS, BOUNDARY_BOTTOM

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		STRUCTURE : SMAP
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Base.@kwdef mutable struct SMAPSTRUCT # <>=<>=<>=<>=<>=<>=<>=<>=<>
            IsTopsoil              :: Vector{Int64}    = Int64[]
            Ks_Impermeable         :: Vector{Float64}  = Float64[]
            RockClass              :: Vector{String}   = String[]
            RockFragment           :: Vector{Float64}  = Float64[]
            Smap_Depth             :: Vector{Float64}  = Float64[]
            Smap_ImpermClass       :: Vector{String}   = String[]
            Smap_MaxRootingDepth   :: Vector{Float64}  = Float64[]
            Smap_PermeabilityClass :: Vector{String}   = String[]
            Smap_SmapFH            :: Vector{String}   = String[]
            Soilname               :: Vector{String}   = String[]
            OrganicMatter          :: Vector{Float64}  = Float64[]
			end
		# ............................................................


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION: READING SMAP                    
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function SMAP(IdSelect, NiZ, path; RockFragment_Max=0.9)
				println("    ~  $(path.inputSmap.Smap) ~")

				Data = CSV.read(path.inputSmap.Smap, DataFrame, header=true)
				DataFrames.sort!(Data, [:Id])

				# reading the structure defined above
				smap = SMAPSTRUCT()
				
            smap.IsTopsoil              = convert(Vector{Float64}, Data."IsTopsoil")
            smap.IsTopsoil              = Int64.(smap.IsTopsoil)
            smap.IsTopsoil              = smap.IsTopsoil[IdSelect]

            smap.RockClass              = convert(Vector{String}, Data."RockClassRock")
            smap.RockClass              = smap.RockClass[IdSelect]
				
            smap.RockFragment           = convert(Vector{Float64}, Data."Stone_Prop")
            smap.RockFragment           = min.(RockFragment_Max, smap.RockFragment)
            smap.RockFragment           = smap.RockFragment[IdSelect]
				
            smap.Smap_Depth             = convert(Vector{Float64}, Data."depth_mm")
            smap.Smap_Depth             = smap.Smap_Depth[IdSelect]
            
            smap.Smap_ImpermClass       = convert(Vector{String}, Data."FHImpermClass")
            smap.Smap_ImpermClass       = smap.Smap_ImpermClass[IdSelect]
				
            smap.Smap_MaxRootingDepth   = convert(Vector{Float64}, Data."MaxRootingDepth_mm")
            smap.Smap_MaxRootingDepth   = smap.Smap_MaxRootingDepth[IdSelect]
				
            smap.Smap_PermeabilityClass = convert(Vector{String}, Data."SoilImpermClass")
            smap.Smap_PermeabilityClass = smap.Smap_PermeabilityClass[IdSelect]

            smap.Smap_SmapFH            = convert(Vector{String}, Data."SmapFH")
            smap.Smap_SmapFH            = smap.Smap_SmapFH[IdSelect]

            smap.Soilname               = convert(Vector{String}, Data."Soilname")
            smap.Soilname               = smap.Soilname[IdSelect]

            smap.OrganicMatter          = convert(Vector{Float64}, Data."OrganicMatter")
            smap.OrganicMatter          = smap.OrganicMatter[IdSelect]

				NiZ = length(smap.Smap_ImpermClass)

				# CORRECTION FOR RockFragment if no correction is required
					for iZ=1:NiZ
						if smap.RockClass[iZ] == "NoAdj"
							smap.RockFragment[iZ] == 0.0
						end
					end

				# KS_IMPERMEABLE
					# If impermeable than the value is greater than 0
					KsImpClass_Dict = readSmap.IMPERMEABLE_CLASS(path.inputSmap.LookupTable_Impermeable)

					smap.Ks_Impermeable = fill(-1.0, NiZ)
					for iZ=1:NiZ
						smap.Ks_Impermeable[iZ] = KsImpClass_Dict[smap.Smap_ImpermClass[iZ]]
					end
			
			return smap
			end  # function: SMAP

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : BOUNDARY_BOTTOM
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function BOUNDARY_BOTTOM(path)

				# Read data
					Data           = CSV.read(path.inputSmap.SoilProfile, DataFrame, header=true)
               BoundarySmap   = convert(Vector{String}, Data."Boundary")
					
               Soilname       = convert(Vector{String}, Data."Soilname")
					
               N_BoundarySmap = length(BoundarySmap)

				# Dictionary of boundary
					Dict_Boundary_Smap2Hypix = Dict{String, String}()

					Terminology_SmapBoundary = ("FreeDrainage","LowImpRock","PermSoil_Fluid","HighImpRock")
					Terminology_HypixBoundary = ("Free","Impermeable","WaterTable","Impermeable")

					for i=1:length(Terminology_SmapBoundary)
						Dict_Boundary_Smap2Hypix[Terminology_SmapBoundary[i]] = Terminology_HypixBoundary[i]
					end

				# Dictionary of names-> Hypix
					Dict_SoilNames_2_HypixBottomBoundary = Dict{String, String}()

				# Converting Smap boundary -> Hypix 
					Hypix_BottomBoundary = fill("", N_BoundarySmap)
					for i = 1:N_BoundarySmap
					# println(BoundarySmap[i])
						Hypix_BottomBoundary[i] = Dict_Boundary_Smap2Hypix[BoundarySmap[i]]

						Dict_SoilNames_2_HypixBottomBoundary[Soilname[i]] = Hypix_BottomBoundary[i] 
					end
				
			return Dict_SoilNames_2_HypixBottomBoundary
			end  # function: BOUNDARY_BOTTOM
		# ------------------------------------------------------------------ 


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
               Data             = CSV.read(Path,  DataFrame, header=true)

               RockClass        = convert(Vector{String}, Data. "RockClass")

               N_RockClass      = length(RockClass)

               RockClass_Unique = unique(RockClass)
					
               N_RockClass      = length(RockClass_Unique)

				# Dictionary
					RockClass_Dict = Dict("a"=>9999)
					for i=1:N_RockClass
						RockClass_Dict[RockClass_Unique[i]] = i
					end

				# Read data
               Ψ₂   = convert(Vector{Float64}, Data."H[mm]")
               θ₂   = convert(Vector{Float64}, Data."Theta[0-1]")

               N₂   = length(Ψ₂)

               Ψ_Rf = zeros(Int64,(N_RockClass, 100))
               θ_Rf = zeros(Float64,(N_RockClass, 100))
               N_Ψ  = zeros(Int64,(N_RockClass))

					iRockClass=1 ; iΨ=1
					for i=1:N₂
						if RockClass[i] == RockClass_Unique[iRockClass]
							Ψ_Rf[iRockClass,iΨ] = Ψ₂[i]
							θ_Rf[iRockClass,iΨ] = θ₂[i]
						else
							N_Ψ[iRockClass]  = iΨ - 1
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
					Data = CSV.read(Path, DataFrame, header=true)

               ImpermClass    = convert(Vector{String}, Data."FHImpermClass")
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



		# ------------------------------------------------------------------
		
	end  # module: smap
	# ...........................................................