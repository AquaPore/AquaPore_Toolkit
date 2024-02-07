# =============================================================
#		module: tableSmap
# =============================================================
module tableSmap
   import ..tool, ..wrc, ..kunsat, ..cst
   import CSV, Tables

   	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         function θΨK(hydro, hydroOther, IdSelect, KₛModel, NiZ, Path, smap)
            println("    ~  $(Path) ~")

            Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(NiZ, hydro)

            Matrix2, FieldName_String2 = tool.readWrite.STRUCT_2_FIELDNAME(NiZ, hydroOther)

            # Concatenating matrices
               Matrix = hcat(Matrix, Matrix2)

               Matrix = hcat(Matrix, KₛModel)

               FieldName_String = vcat(FieldName_String, FieldName_String2)

               FieldName_String = vcat("Id", "SoilName", "Depth", FieldName_String, "Ks_Model[mm/s]")

               CSV.write(Path, Tables.table( [string.(IdSelect[1:NiZ]) smap.Soilname[1:NiZ] smap.Smap_Depth[1:NiZ] Matrix]), writeheader=true, header=FieldName_String, bom=true)
         return nothing
         end  # function:  θΨK

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : Smap
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   """
   TopnetModel = ["ThetaS_Ch[mm3_mm3]";"ThetaR_Ch[mm3_mm3]";"LambdaCh_Ch[-]";"Hch_Ch[mm]";" Hga_Ch[mm]";"Ks_Vg[mm_s1]; "0mm"; "500mm"; "1000mm"; "2000mm"; "4000mm"; "10000mm"; "150000mm"]

   JulesModel_CH =  ["ThetaS_Ch[mm3_mm3]";"ThetaR_Ch[mm3_mm3]";"LambdaCh_Ch[-]";"Hch_Ch[mm]"; "Ks_Ch[mm_s1]";" Hga_Ch[mm]";"3300mm";"10000mm" ;"150000mm"]

   JulesModel_VangenuchtenJules = ["ThetaS_VgJules[mm3_mm3]";"ThetaR_VgJules[mm3_mm3]";"n_VgJules[-]";"Hvg_VgJules[mm]"; "Ks_VgJules[mm_s1]";"3300mm";"10000mm"]

   """
      function SMAP(hydro, IdSelect, IsTopsoil, NiZ, optionₘ, param, path, smap)

         println("    ~  $(path.tableSmap.Table_Smap) ~")

          # User input
            HeaderSmap = true # <true> the greek characters are replaced by alphabet; <false> original parameter names with no units usefull to use values in SoilWater-ToolBox

            🎏_BrooksCorey       = true
            🎏_ClappHornberger   = true
            🎏_VanGenuchten      = false
            🎏_VanGenuchtenJules = false
            🎏_Kosugi            = true
            🎏_Kosugi_Table_θψ   = true
            🎏_Kosugi_Table_Kψ   = true
            🎏_Fc_Pwp_Paw        = true

         Header = ["Id"; "SoilName"; "Depth_mm"; "IsTopsoil"; "RockFragment_%";"MaxRootingDepth_mm"; "PermeabilityClass"; "SmapFH"; "ImpermClass"]
         Data = []
      
      # Select data
         # HydroModel_θΨ == "BrooksCorey"
         if 🎏_BrooksCorey # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "BrooksCorey"

            Path_θΨ =  path.tableSoilwater.Path_Soilwater_Table *  "_" * string(HydroModel_θΨ) *  "_" * "Table_Smap_θΨK.csv"
            
            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"λbc";"Ψbc"; "Ks"; "Ψga"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))

               try
                  Data = hcat(Data[1:NiZ, :], Data_θΨ[1:NiZ, :])
               catch
                  Data = Data_θΨ[1:NiZ, :]
               end
               
               if HeaderSmap
                  Header_θΨ = ["ThetaS_Bc[mm3_mm3]";"ThetaR_Bc[mm3_mm3]";"LambdaBc_Bc[-]";"Hbc_Bc[mm]";"Ks_Bc[mm_s1]";"Hga_Bc[mm]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # 🎏_BrooksCorey


         # HydroModel_θΨ == "ClappHornberger"
         if  🎏_ClappHornberger # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "ClappHornberger"

            Path_θΨ =  path.tableSoilwater.Path_Soilwater_Table *  "_" * string(HydroModel_θΨ) *  "_" * "Table_Smap_θΨK.csv"

            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"λch";"Ψch";"Ks";"Ψga"]

               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))

               try
                  Data = hcat(Data[1:NiZ, :], Data_θΨ[1:NiZ, :])
               catch
                  Data = Data_θΨ[1:NiZ, :]
               end

               if HeaderSmap
                  Header_θΨ = ["ThetaS_Ch[mm3_mm3]";"ThetaR_Ch[mm3_mm3]";"LambdaCh_Ch[-]";"Hch_Ch[mm]"; "Ks_Ch[mm_s1]";" Hga_Ch[mm]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end #  🎏_ClappHornberger

         
         # HydroModel_θΨ == "Vangenuchten"
         if 🎏_VanGenuchten # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "Vangenuchten"

            Path_θΨ =  path.tableSoilwater.Path_Soilwater_Table *  "_" * string(HydroModel_θΨ) *  "_" * "Table_Smap_θΨK.csv"

            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"N";"Ψvg"; "Ks"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:NiZ, :], Data_θΨ[1:NiZ, :])
               catch
                  Data = Data_θΨ[1:NiZ, :]
               end
         
               if HeaderSmap
                  Header_θΨ = ["ThetaS_Vg[mm3_mm3]";"ThetaR_Vg[mm3_mm3]";"N_Vg[-]";"Hvg_Vg[mm]"; "Ks_Vg[mm_s1]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # 🎏_VanGenuchten


         # HydroModel_θΨ == "VangenuchtenJules"
         if 🎏_VanGenuchtenJules # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "VangenuchtenJules"

            Path_θΨ =  path.tableSoilwater.Path_Soilwater_Table *  "_" * string(HydroModel_θΨ) *  "_" * "Table_Smap_θΨK.csv"

            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"N";"Ψvg"; "Ks"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:NiZ, :], Data_θΨ[1:NiZ, :])
               catch
                  Data = Data_θΨ[1:NiZ, :]
               end
         
               if HeaderSmap
                  Header_θΨ = ["ThetaS_VgJules[mm3_mm3]";"ThetaR_VgJules[mm3_mm3]";"n_VgJules[-]";"Hvg_VgJules[mm]"; "Ks_VgJules[mm_s1]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # 🎏_VanGenuchten


         # HydroModel_θΨ == "Kosugi"
         if 🎏_Kosugi # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "Kosugi"

            Path_θΨ =  path.tableSoilwater.Path_Soilwater_Table *  "_" * string(HydroModel_θΨ) *  "_" * "Table_Smap_θΨK.csv"

            if isfile(Path_θΨ)
               Select_θΨ =["θs";"θr";"Ks";"Ψm";"σ";"θsMacMat_ƞ";"σMac";"ΨmMac"; "θsMacMat";"Φ"]
                   
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:NiZ, :], Data_θΨ[1:NiZ, :])
               catch
                  Data = Data_θΨ[1:NiZ, :]
               end
         
               if HeaderSmap
                  Header_θΨ = ["ThetaS_Kg[mm3_mm3]";"ThetaR_Kg[mm3_mm3]";"Ks_Kg[mm_s1]";"Hm_Kg[mm]";"Sigma_Kg";"ThetaSMacMatNorm_Kg[mm3_mm3]";"SigmaMac_Kg";"HmMac_Kg[mm]";"ThetaSMacMat_Kg[mm3_mm3]";"TotalPorosity[mm3_mm3]"]
               else
                  Header_θΨ = Select_θΨ
               end
               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)

         end # 🎏_Kosugi


        

      # CREATING TABLES θ(ψ) & TABLES K(ψ)==========
         Path_Select_θΨ = path.tableSoilwater.Path_Soilwater_Table *  "_Kosugi_Table_Smap_Select_θΨ.csv"
         Path_Select_KΨ = path.tableSoilwater.Path_Soilwater_Table *  "_Kosugi_Table_Smap_Select_KΨ.csv"
         
         if optionₘ.HydroModel⍰ == "Kosugi"
            # Creating Table θ(Ψ)
               N_Ψ = length(param.smap.Ψ_Table[:])
               θ₂  = fill(0.0::Float64, (NiZ, N_Ψ))

               for iZ=1:NiZ
                  for iΨ =1:N_Ψ
                     Ψ₂         = param.smap.Ψ_Table[iΨ]
                     θ₂[iZ, iΨ] = wrc.Ψ_2_θ(optionₘ, Ψ₂, iZ, hydro)
                  end # iΨ
               end # iZ
            
               Header_θΨ = "ThetaH_" .* string.(Int64.(param.smap.Ψ_Table)) .* "_mm"  
               CSV.write(Path_Select_θΨ, Tables.table( θ₂[1:NiZ, :]), writeheader=true, header=Header_θΨ, bom=!(HeaderSmap))

            # Creating Table K(Ψ) ======
               N_Ψ = length(param.hydro.K_Table[:])
               K₂  = fill(0.0::Float64, (NiZ, N_Ψ))

               for iZ=1:NiZ
                  for iΨ =1:N_Ψ
                     Ψ₂ = param.hydro.K_Table[iΨ]
                     K₂[iZ, iΨ] = kunsat.KUNSAT_θΨSe(optionₘ, Ψ₂, iZ, hydro)
                  end # iΨ
               end # iZ

               Header_KΨ = "KunsatH_" .* string.(Int64.(param.hydro.K_Table)) .* "_mm s1" 

               CSV.write(Path_Select_KΨ, Tables.table( K₂[1:NiZ, :]), writeheader=true, header=Header_KΨ, bom=!(HeaderSmap))
         end # CREATING TABLES θ(ψ) & TABLES K(ψ)


      # FIELD CAPACITY, PERMENANT WILTING POINT, AVAILABLE WATER CONTENT 


      # WRITTING DATA TO Table_Smap <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
         if (🎏_Kosugi_Table_θψ || 🎏_Fc_Pwp_Paw) && isfile(Path_Select_θΨ)
               Header_θΨ = "ThetaH_" .* string.(Int64.(param.smap.Ψ_Table .* cst.Mm_2_kPa)) .* "_kPa" 

               θ₂ = Tables.matrix(CSV.File(Path_Select_θΨ))

               Data = [Data[1:NiZ, :] θ₂[1:NiZ, :]]

               Header =  append!(Header, Header_θΨ)

					if 🎏_Fc_Pwp_Paw
                  θ_0kPa        = zeros(NiZ)
                  θ_5kPa        = zeros(NiZ)
                  θ_10kPa       = zeros(NiZ)
                  θ_1500kPa     = zeros(NiZ)
                  Paw           = zeros(NiZ)
                  Macroporosity = zeros(NiZ)
						AirfilledMacroporosity = zeros(NiZ)

						for iZ=1:NiZ
                     iθ_0kPa           = findfirst(isequal(0.0), param.smap.Ψ_Table)
                     θ_0kPa[iZ]        = θ₂[iZ, iθ_0kPa]

							iθ_5kPa           = findfirst(isequal(500.0), param.smap.Ψ_Table)
                     θ_5kPa[iZ]        = θ₂[iZ, iθ_5kPa]

                     iθ_10kPa          = findfirst(isequal(1000.0), param.smap.Ψ_Table)
                     θ_10kPa[iZ]       = θ₂[iZ, iθ_10kPa]

                     iθ_1500kPa        = findfirst(isequal(150000.0), param.smap.Ψ_Table)
                     θ_1500kPa[iZ]     = θ₂[iZ, iθ_1500kPa]

                     Paw[iZ]                = θ_10kPa[iZ] - θ_1500kPa[iZ]
                     Macroporosity[iZ]      = θ_0kPa[iZ] - θ_5kPa[iZ]
                     AirfilledMacroporosity[iZ] = θ_0kPa[iZ] - θ_10kPa[iZ]
						end

						Header_Pc_Pwp_Paw = ["Fc", "Pwp", "Paw", "Macroporosity", "AirfilledMacroporosity"]
						Header =  append!(Header, Header_Pc_Pwp_Paw)

						Data = [Data[1:NiZ, :] θ_10kPa[1:NiZ] θ_1500kPa[1:NiZ] Paw[1:NiZ] Macroporosity[1:NiZ] AirfilledMacroporosity[1:NiZ]]
					end
         end # 🎏_Kosugi

         if 🎏_Kosugi_Table_Kψ  && isfile(Path_Select_KΨ)
            K₂= Tables.matrix(CSV.File(Path_Select_KΨ))

            Data = [Data[1:NiZ, :] K₂[1:NiZ, :]]
      
            Header_KΨ = "KunsatH_" .* string.(round.(param.hydro.K_Table .* cst.Mm_2_kPa, digits=2)) .* "kPa_[mm s-1]"            

            Header =  append!(Header, Header_KΨ)
         end # 🎏_Kosugi


      # COMBINING OUTPUTS =====================================================  
         CSV.write(path.tableSmap.Table_Smap, Tables.table( [IdSelect[1:NiZ] smap.Soilname[1:NiZ] smap.Smap_Depth[1:NiZ] smap.IsTopsoil[1:NiZ] smap.RockFragment[1:NiZ] smap.Smap_MaxRootingDepth[1:NiZ] smap.Smap_PermeabilityClass[1:NiZ] smap.Smap_SmapFH[1:NiZ] smap.Smap_ImpermClass[1:NiZ] Data[1:NiZ,:]]), writeheader=true, header=Header, bom=!(HeaderSmap))

      return nothing
      end  # function:  smap
	# ............................................................
  
end  # module: tableSmap
# ............................................................