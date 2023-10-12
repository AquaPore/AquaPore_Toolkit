# =============================================================
#		module: plotSmap
# =============================================================
module plotSmap

      # =============================================================
      #		module: makie
      # =============================================================
      module makie
         import ...cst, ...hydroStruct, ...kunsat, ...reading, ...readSmap, ...wrc, ...tool
         using CairoMakie
      
         # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         #		FUNCTION : HYDROPARAM
         # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         function HYDROPARAM(hydro, IdSelect, K_KΨobs, KₛModel, NiZ, N_KΨobs, N_θΨobs, option, path, Smap_Depth, Soilname, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs; N_Se=1000)
            println("  ==  START: Plotting HydroParam  ==")
      
            🎏_OtherData1 = false
            🎏_OtherData2 = false
      
            # ===================== DATA =====================
            θ_Sim             = fill(0.0,N_Se)
            θ_OtherData       = fill(0.0,N_Se)
            θ_OtherData2      = fill(0.0,N_Se)
            θ_OtherData3      = fill(0.0,N_Se)
            Kunsat_Sim        = fill(0.0,N_Se)
            Kunsat_OtherData  = fill(0.0,N_Se)
            Kunsat_OtherData2 = fill(0.0,N_Se)
            Kunsat_OtherData3 = fill(0.0,N_Se)
            θobs =[]

            Ψ_θΨobs_Min = 0.0

            if 🎏_OtherData1
               Path = "D:\\Main\\MODELS\\SoilWater-ToolBox2\\src\\INPUT\\Data_SoilWater\\Smap20210226\\Smap20210226_ClappHornberger_Constrained_A_Table_ThetaHK.csv"

               option.hydro.HydroModel⍰ = "ClappHornberger"
                  # Structure of the hydroparameters
                     hydroData = hydroStruct.HYDROSTRUCT(option.hydro, NiZ)  
                  # Populate the values of the parameters
                  option.hydro.HydroModel⍰ = "ClappHornberger"
                     hydroData, ~ = tool.readWrite.READ_STRUCT_SIMPLE(hydroData, Path)

               Path = "D:\\Main\\MODELS\\SoilWater-ToolBox2\\src\\INPUT\\Data_SoilWater\\Smap20210226\\Smap20210226_Loam.csv"
                  option.hydro.HydroModel⍰ = "ClappHornberger"
                  # Structure of the hydroparameters
                     hydroData2 = hydroStruct.HYDROSTRUCT(option.hydro, NiZ)
                  # Populate the values of the parameters
                     hydroData2, ~ = tool.readWrite.READ_STRUCT_SIMPLE(hydroData2, Path) 

                 Path =  "D:\\Main\\MODELS\\SoilWater-ToolBox2\\src\\INPUT\\Data_SoilWater\\Smap20210226\\Smap20210226_VangenuchtenJules_Constrained_A_Table_ThetaHK.csv"
                  option.hydro.HydroModel⍰ = "VangenuchtenJules"
                  # Structure of the hydroparameters
                     hydroData3 = hydroStruct.HYDROSTRUCT(option.hydro, NiZ)
                  # Populate the values of the parameters
                     hydroData3, ~ = tool.readWrite.READ_STRUCT_SIMPLE(hydroData3, Path) 

            end # if 🎏_OtherData

            for iZ = param.globalparam.N_iZ_Plot_Start:param.globalparam.N_iZ_Plot_End
               Ψ_θΨobs_Max = maximum(Ψ_θΨobs[iZ,N_θΨobs[iZ]]) + 100000.0

               Ψ_Sim = expm1.(range(log1p(Ψ_θΨobs_Min), stop=log1p(Ψ_θΨobs_Max), length=N_Se)) 

               θ_θΨobs_Max = hydro.Φ[iZ]

               # Simulated 
                  for iΨ = 1:N_Se
                     option.hydro.HydroModel⍰ = "Vangenuchten"
                      θ_Sim[iΨ] = wrc.Ψ_2_θ(option.hydro,Ψ_Sim[iΨ], iZ, hydro)
                      Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(option.hydro, Ψ_Sim[iΨ], iZ, hydro)

                     if 🎏_OtherData1
                        # ClappHornberger model Smap_Hydro
                        option.hydro.HydroModel⍰ = "ClappHornberger"
                           θ_OtherData[iΨ] = wrc.Ψ_2_θ(option.hydro,Ψ_Sim[iΨ], iZ, hydroData)
                           option.hydro.HydroModel⍰ = "ClappHornberger"
                           Kunsat_OtherData[iΨ] = kunsat.Ψ_2_KUNSAT(option.hydro, Ψ_Sim[iΨ], iZ, hydroData)

                        # ClappHornberger Loam
                        option.hydro.HydroModel⍰ = "ClappHornberger"
                           θ_OtherData2[iΨ] = wrc.Ψ_2_θ(option.hydro,Ψ_Sim[iΨ], iZ, hydroData2)
                            option.hydro.HydroModel⍰ = "ClappHornberger"
                           Kunsat_OtherData2[iΨ] = kunsat.Ψ_2_KUNSAT(option.hydro, Ψ_Sim[iΨ], iZ, hydroData2)
          
                        # VanGenuchten_Jules
                        option.hydro.HydroModel⍰ = "VangenuchtenJules"
                         θ_OtherData3[iΨ] = wrc.Ψ_2_θ(option.hydro,Ψ_Sim[iΨ], iZ, hydroData3)
                         option.hydro.HydroModel⍰ = "VangenuchtenJules"
                        Kunsat_OtherData3[iΨ] =  kunsat.Ψ_2_KUNSAT(option.hydro, Ψ_Sim[iΨ], iZ, hydroData3)

                         θobs₀ =[ [ 0.456,	0.35,	0.28,	0.16],
                        [0.4465,	0.32,	0.25,	0.15],
                        [0.4465,	0.3,	0.17,	0.1],
                        [0.646,	0.46,	0.36,	0.25],
                        [0.6745,	0.48,	0.36,	0.27],
                        [0.6935,	0.52,	0.43,	0.29]]

                         θobs = θobs₀[iZ,:][1] 
                     else
                         option.hydro.HydroModel⍰ = "VangenuchtenJules"
                         θ_Sim[iΨ] = wrc.Ψ_2_θ(option.hydro,Ψ_Sim[iΨ], iZ, hydro)

                         option.hydro.HydroModel⍰ = "VangenuchtenJules"
                        Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(option.hydro, Ψ_Sim[iΨ], iZ, hydro)
                     end
                  end # iΨ

               # == Title == 
                  Title = Soilname[iZ]  * "_" * string(Int64(floor(Smap_Depth[iZ]))) * "_" * string(option.hydro.HydroModel⍰)

                  Title = Soilname[iZ]  * "  " * string(Int64(floor(Smap_Depth[iZ]))) * " mm"
                  # Title = Title  * "_" * string(option.hydro.σ_2_Ψm⍰)

      
               #  ===================== PLOTTING =====================
                  Fig = Figure()
                              
               #  == Plot_θ_Ψ  ==
                  # Plot_θ_Ψ: General attributes

                  Fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),
                     resolution = (1000, 700))

                  Axis1 = Axis(Fig[1,1], resolution = (1000, 700))

                  xlims!(Axis1, log1p.(cst.Mm_2_kPa * Ψ_θΨobs_Min), log1p.(cst.Mm_2_kPa * Ψ_θΨobs_Max * 1.1))
                  ylims!(Axis1, 0.0, 0.75)
   
                  Axis1.xticks = (log1p.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]]), string.(Int64.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]])))
                  Axis1.xlabel = "ln(1 + Ψ) [kPa]"
                  Axis1.ylabel =  "θ [mm³ mm⁻³]"
                  Axis1.title = Title
                  Axis1.titlesize= 24

                  P_Smap= scatter!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_θΨobs[iZ,1:N_θΨobs[iZ]]), Float64.(θ_θΨobs[iZ,1:N_θΨobs[iZ]]), color=:blue, markersize=15, marker = '■', label="Smap")

                  # Plot_θ_Ψ: Simulated
                     # lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_Sim[1:N_Se], color=:blue, linewidth=2, label="Sim")

                  if 🎏_OtherData1
                     P_ClappHonb = lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_OtherData[1:N_Se], color=:red, linewidth=3)

                     P_ClappHonb_Loan = lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_OtherData2[1:N_Se], color=:yellow1, linewidth=3)

                     P_vangJules = lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_OtherData3[1:N_Se], color=:green, linewidth=3)

                     Ψ = [0.0 ,1000.0,10000.0,150000.0]
                     P_Lab = scatter!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ[1:4]), θobs[1:4], color=:darkviolet, linewidth=2, markersize=15)

                     P_Vang = lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_Sim[1:N_Se], color=:pink, linewidth=3)
                  else
                      lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_Sim[1:N_Se], color=:blue, linewidth=2)

                  end

                  # Plot_θ_Ψ: Total porosity point
                     X = zeros(Float64,1)
                     X[1] = 0.0
                     Y = zeros(Float64,1)
                     Y[1] = hydro.Φ[iZ]
                     Label = 
                     P_PtotalPorosity = scatter!(Fig[1,1], log1p.(cst.Mm_2_kPa .* X), Y, color=:slateblue3, markersize=20, marker ="●")

               # == Plot_K_Ψ  ==
               if "Ks" ∈ optim.ParamOpt
                     Axis2 = Axis(Fig[1,2])
                     Axis2.xticks = (log1p.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]]), string.(Int64.(cst.Mm_2_kPa * Ψ_θΨobs[iZ,1:N_θΨobs[iZ]])))
                     Yticks = 1:1:6
                     Axis2.yticks = (Yticks,string.(Yticks))
                     Axis2.xlabel = "ln(1 + Ψ) [kPa]"
                     Axis2.ylabel =  "ln ( 1 + K (Ψ) ) [mm h⁻¹]"

                   Pvang = lines!(Fig[1,2], log1p.(Ψ_Sim[1:N_Se].*cst.Mm_2_kPa), log1p.(Kunsat_Sim[1:N_Se].*cst.MmS_2_MmH), color=:pink, linewidth=3)

                  if 🎏_OtherData1

                     Pclapp = lines!(Fig[1,2], log1p.(Ψ_Sim[1:N_Se].*cst.Mm_2_kPa), log1p.(Kunsat_OtherData[1:N_Se].*cst.MmS_2_MmH), color=:red, linewidth=3)

                     Ploan = lines!(Fig[1,2], log1p.(Ψ_Sim[1:N_Se].*cst.Mm_2_kPa), log1p.(Kunsat_OtherData2[1:N_Se].*cst.MmS_2_MmH), color=:yellow1, linewidth=3)

                     P9 = lines!(Fig[1,2], log1p.(Ψ_Sim[1:N_Se].*cst.Mm_2_kPa), log1p.(Kunsat_OtherData3[1:N_Se].*cst.MmS_2_MmH), color=:green, linewidth=3)

                  else
                      Axis2 = Axis(Fig[1,2])

                     K_Ψ_Max = maximum(K_KΨobs[iZ,1:N_KΨobs[iZ]])
                     xlims!(Axis2, log1p.(cst.Mm_2_kPa*Ψ_θΨobs_Min), log1p.(Ψ_θΨobs_Max*cst.Mm_2_kPa))
                     ylims!(Axis2,  (log1p(0.0), log1p(K_Ψ_Max* cst.MmS_2_CmH * 1.1)))
                     Axis2.xlabel = "ln(1 + Ψ) [kPa]"
                     Axis2.ylabel =  "ln ( 1 + K (Ψ) ) [cm h⁻¹]"
                       X = zeros(Float64,1)
                     X[1] = 0.0 
                     Y = zeros(Float64,1)
                     Y[1] = K_Ψ_Max
                     Label = "Ks_Max"
                     scatter!(Fig[1,2], log1p.(X.* cst.Mm_2_kPa) , log1p.(Y.*cst.MmS_2_CmH) ,Y, color=:yellow, markersize=15, marker = '■', label=Label )

                  # PlotK_Ψ: K(Ψ) obs
                     X = Ψ_KΨobs[iZ,1:N_KΨobs[iZ]]
                     Y = K_KΨobs[iZ,1:N_KΨobs[iZ]]
                     Label = "Obs"
                     scatter!(Fig[1,2], log1p.(X.* cst.Mm_2_kPa) , log1p.(Y.*cst.MmS_2_CmH) ,Y, color=:red, markersize=10, marker = '■', label=Label )

                  # Plot_K_Ψ: K(Ψ) sim
                     X = Ψ_Sim
                     Y = Kunsat_Sim 
                     Label = "Sim"
                     lines!(Fig[1,2], log1p.(Ψ_Sim.*cst.Mm_2_kPa), log1p.(Kunsat_Sim.*cst.MmS_2_CmH), color=:blue, label=Label)
                     
                  end
                  
               end # "Ks" ∈ optim.ParamOpt 

               # Fig[3, 1] = Legend(Fig, Axis1, "PLOTS", orientation=:horizontal)

               # Path = path.plotSoilwater.Plot_θΨK * "Lab_ThetaH_" * Title * ".svg" 

               leg = Fig[1, end+1] = Legend(Fig, [P_Smap, P_ClappHonb, P_ClappHonb_Loan, P_vangJules, P_Vang, P_Lab, P_PtotalPorosity], ["Smap", "ClapHornberger", "ClapHornberger_Loam", "VangJules", "Vang","Lab","TotalPorosity"])

               Fig[2, 1:2] = leg
               trim!(Fig.layout)
               leg.orientation = :horizontal
               trim!(Fig.layout)
               leg.tellheight = true
               
               Path = path.plotSoilwater.Plot_θΨK * "Lab_ThetaH_" * string(path.option.ModelName) * "_" * string(IdSelect[iZ]) * ".svg" 
      
               save(Path, Fig)
     
               # Displaying figure in VScode
               if option.general.PlotVscode
                  display(Fig)
               end

               println("    ~  $(Path) ~")
            end # for iZ

            println("  ==  END: Plotting HydroParam  == \n")		
         return nothing
         end  # function: HYDROPARAM
            
      end  # module: makie
      # ............................................................
   
end  # module: plotSmap
# ............................................................