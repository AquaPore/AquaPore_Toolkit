# =============================================================
#		module: smap2hypix
# =============================================================
module smap2hypix
   import  ..hydroStruct, ..tool, ..wrc
   using Tables, CSV, DataFrames

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : SMAP_2_HYDRO
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function SMAP_2_HYPIX(Dict_SoilNames_2_HypixBottomBoundary, hydro, NiZ, optionₘ, param, path, RockFragment, smap)
    
         # Index of each soil profile in the spreadsheet
            iSoilProfile_End, iSoilProfile_Start, N_SoilProfile, Soilname_SoilProfile = SMAP_SOILPROFILE(NiZ, smap)

         # PATHS
            Path_OutputHydro = joinpath(path.smap2Hypix.Path_Smap2Hypix, "HYDRO_INPUT")
               if isfile(Path_OutputHydro)
                  rm(Path_OutputHydro; force=true, recursive=true)
               else
                  mkpath(Path_OutputHydro)
               end

            Path_OutputSoilHorizon = joinpath(path.smap2Hypix.Path_Smap2Hypix, "SOILLAYER")
               if isfile(Path_OutputSoilHorizon)
                  rm(Path_OutputSoilHorizon; force=true, recursive=true)
               else
                  mkpath(Path_OutputSoilHorizon)
               end

            Path_OutputExtra = joinpath(path.smap2Hypix.Path_Smap2Hypix, "SOILPROFILE")
               if isfile(Path_OutputExtra)
                  rm(Path_OutputExtra; force=true, recursive=true)
               else
                  mkpath(Path_OutputExtra)
               end

         # For every soil profile
         for iSoilProfile=1:N_SoilProfile
            # iStart and iEnd of Horizon
               iHorizon_Start= iSoilProfile_Start[iSoilProfile]
               iHorizon_End = iSoilProfile_End[iSoilProfile]

            # Checking
               iHorizon_End= min(NiZ, iHorizon_End)
            
					N_Horizon = iHorizon_End -  iSoilProfile_Start[iSoilProfile] + 1 # Final

            # WRITING HYDRAUliC PARAM TO HYPIX
               hydroSmap = hydroStruct.HYDROSTRUCT(optionₘ, N_Horizon) # Making a structure

               # Putting the values of Output_Vector into structure
               for iFieldname in propertynames(hydro)
                  Soil_Horizon = getfield(hydro, iFieldname)[iHorizon_Start:iHorizon_End]		
                  setfield!(hydroSmap, Symbol(iFieldname), Soil_Horizon)
               end

             # Writing to file
               Path_OutputHydro₂ = joinpath(Path_OutputHydro, Soilname_SoilProfile[iSoilProfile] *  "_TableHydro.csv")
               TABLE_HYDRO(hydroSmap, N_Horizon, Path_OutputHydro₂)
            
            # WRITING DEPTH OF SOIL
               Zhorizon = smap.Smap_Depth[iHorizon_Start:iHorizon_End]

               Path_OutputSoilHorizon₂ = joinpath(Path_OutputSoilHorizon, Soilname_SoilProfile[iSoilProfile] *  "_SoilLayer.csv") 
               Zhorizon = LAYER_DISCRETISATION(hydroSmap, N_Horizon, Path_OutputSoilHorizon₂, Zhorizon)

            # WRITTING MAXIMUM ROOTING DEPTH
             #  The maximum rooting depth is a repeat between iHorizon_Start:iHorizon_End
               RootingDepth = min(smap.Smap_MaxRootingDepth[iHorizon_Start], maximum(Zhorizon))

            # Maximum Rock Fragment of the soil profile
               # MAximum RockFragment for the top 400mm
              Nz_400 = iHorizon_End
               for iZ = iHorizon_Start:iHorizon_End
                  if smap.Smap_Depth[iZ] > 400.0
                     Nz_400 = iZ
                     break
                  end
               end # iZ=1:Nz

               Rf_Max = maximum(RockFragment[iHorizon_Start:Nz_400])

            # Boundary conditions
               HypixBottomBoundary = Dict_SoilNames_2_HypixBottomBoundary[Soilname_SoilProfile[iSoilProfile]]

            # WRITTING TO SOILPROFILE FILE
               Path_OutputExtra₂ = joinpath(Path_OutputExtra, Soilname_SoilProfile[iSoilProfile] *  "_Soilprofile.csv")
               SOILPROFILE(HypixBottomBoundary, Path_OutputExtra₂, Rf_Max, RootingDepth)   
         end #  iSoilProfile=1:N_SoilProfile      
   
      return nothing
      end  # function: SMAP_2_HYDRO
   #----------------------------------------------------------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : SMAP_SOIL_iSoilProfile
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function SMAP_SOILPROFILE(NiZ, smap)
         iSoilProfile_End = []
         iSoilProfile_Start = [1]
         Soilname_Initial = smap.Soilname[1]
         Soilname_SoilProfile = [smap.Soilname[1]]
         i = 1
         N_SoilProfile = 1
         for iSoilname in smap.Soilname
            # if soil changes
            if iSoilname ≠ Soilname_Initial
               append!(iSoilProfile_Start, i)
               append!(iSoilProfile_End, i-1)
               push!(Soilname_SoilProfile, iSoilname)

               Soilname_Initial = smap.Soilname[i] # New soil
               N_SoilProfile += 1
            elseif  i == NiZ
               append!(iSoilProfile_End, i)  
            end  # if: name
            i += 1
         end # iSoilname
         
      return iSoilProfile_End, iSoilProfile_Start, N_SoilProfile, Soilname_SoilProfile
      end  # function: SMAP_SOIL_iSoilProfile
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : LAYER_DISCRETISATION
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function LAYER_DISCRETISATION(hydroSmap, N_Horizon, Path_OutputSoilHorizon₂, Zhorizon; Seᵢₙᵢ=0.6)

         # COMPUTING θini
            θᵢₙᵢ = fill(0.0::Float64, N_Horizon)

            # Assuming that all SoilProfiles have the same Se
            for iZ=1:N_Horizon
               θᵢₙᵢ[iZ] = wrc.Se_2_θ(Se₁=Seᵢₙᵢ, θs=hydroSmap.θs[iZ], θr=hydroSmap.θr[iZ])
            end

         # WRITTING TO FILE
            Header = ["Layer";"Z"; "θini"]

            iZ = collect(1:1:length(Zhorizon))
   
            CSV.write(Path_OutputSoilHorizon₂, Tables.table([Int64.(iZ) Zhorizon θᵢₙᵢ]), writeheader=true, header=Header, bom=true)
            
      return Zhorizon
      end  # function: LAYER_DISCRETISATION
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TABLE_HYDRO
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TABLE_HYDRO(hydroSmap, N_iSoilProfiles, Path)
			# println("			~ $(Path) ~")

			Id = 1:1:N_iSoilProfiles

			Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_iSoilProfiles, hydroSmap)
					
			pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

         CSV.write(Path, Tables.table([Int64.(Id) Matrix]), writeheader=true, header=FieldName_String, bom=true)
      return nothing
		end  # function: TABLE_HYDRO
   #...................................................................  


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : SOILPROFILES
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function SOILPROFILE(HypixBottomBoundary::String, Path_OutputExtra₂::String, Rf_Max::Float64, RootingDepth::Float64)
         Header = ["RootingDepth[mm]"; "BottomBoundary"; "Rf_Max[%]"]

         CSV.write(Path_OutputExtra₂, Tables.table([RootingDepth HypixBottomBoundary Rf_Max]), writeheader=true, header=Header, bom=true)
         
      return nothing
      end  # function: EXTRAS
   # ------------------------------------------------------------------

end # module smap2hypix
