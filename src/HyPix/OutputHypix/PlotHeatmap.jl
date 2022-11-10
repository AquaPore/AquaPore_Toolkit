
using DataFrames, CSV, Dates, CairoMakie

function PLOT_HEATMAP()

   θZₐᵥₑᵣ =  [100.0, 200.0, 300.0, 400.0, 500.0, 600.0]

   Option_Data = :YearMonth # :YearMonth :MonthCombine

   SiteName =["Awatere", "Balclutha", "Cromwell", "Dannevirke", "Darfield", "Dargaville", "Dunedin", "Greymouth", "Hamilton", "Hanmer", "LakeTekapo", "Lauder", "Lincoln", "Martinborough", "Middlemarch",  "Paraparaumu", "Ranfurly", "Rangiora", "Reefton", "Stratford",  "Turangi", "Waipawa", "Wallaceville", "Winchmore", "Windsor"]

   for iSiteName in SiteName
      println("<==== ", iSiteName, " ====>")
      PLOT_HEATMAP_2(iSiteName, θZₐᵥₑᵣ, Option_Data)    
   end
end

function PLOT_HEATMAP_2(SiteName, θZₐᵥₑᵣ, Option_Data)
   PathInput₀ = "D:\\Main\\MODELS\\SoilWater_ToolBox\\data\\OUTPUT\\Hypix\\SMAP\\"

   if Option_Data == :YearMonth

      PathInput = PathInput₀ * SiteName * "/Table/" * SiteName * "_Daily_Table_StatisticZaver_YearMonth_1.csv"

      PathOutput = PathInput₀ * SiteName * "/Plots/" * SiteName * "_Heatmap_YearMonth.svg"

   elseif Option_Data == :MonthCombine
      PathInput = PathInput₀ * SiteName * "/Table/" * SiteName * "_Daily_Table_Statistic_MonthCombine_1.csv"

      PathOutput = PathInput₀ * SiteName * "/Plots/" * SiteName * "_Heatmap_MonthCombine.svg"
   else
      error("Option_Data selection not found")
   end

   Df = DataFrame(CSV.File(PathInput))

   if Option_Data == :YearMonth
      Df_θZₐᵥₑᵣ =  "Z" .* string.(θZₐᵥₑᵣ) .*"mm" .* "_mean"
   elseif Option_Data == :MonthCombine
      Df_θZₐᵥₑᵣ =  "Z" .* string.(θZₐᵥₑᵣ) .*"mm" .* "_mean" .* "_mean"
   end

   Header_θZₐᵥₑᵣ = string.(θZₐᵥₑᵣ)

   Header_θZₐᵥₑᵣ = reverse(Header_θZₐᵥₑᵣ)

   Nit = length(Df[!,Df_θZₐᵥₑᵣ[1]])
   Width = max(Nit * 40, 800)

   HEATMAP(Df, Df_θZₐᵥₑᵣ, Header_θZₐᵥₑᵣ, Nit, Option_Data, PathOutput, SiteName; Height= 500, Width = Width)
end

function HEATMAP(Df, Df_θZₐᵥₑᵣ, Header_θZₐᵥₑᵣ, Nit, Option_Data, PathOutput, SiteName; Height=400, Width=1000)

   N_θZₐᵥₑᵣ = length(Df_θZₐᵥₑᵣ)

   if Option_Data == :YearMonth
      Xticks = string.(Df[!,:Year]) .* "_" .* monthabbr.(Df[!,:Month])
   elseif Option_Data == :MonthCombine
      Xticks = string.(Df[!,:MonthlyCombine])
   end

   Title = "©HyPix,    SiteName = $SiteName" 
   Fig = Figure(fontsize=20, backgroundcolor = RGBf(0.98, 0.98, 0.98), font="CMU Serif")

   Ax1 = Axis(Fig[1,1], title=Title, titlesize=40, xticks =(1:Nit, Xticks), ylabel= L"$\Delta Flux$ [mm]", height=Height, width=Width, xgridvisible = true, ygridvisible = false, yticks = LinearTicks(10))
      xlims!(Ax1, 1, Nit)

      hidexdecorations!(Ax1, grid=false, ticks=true, ticklabels=true)
      
      if Option_Data == :YearMonth
         barplot!(Ax1, 1:Nit,Df[!,:ΔPrThroughfall_sum], color=:blue, strokecolor=:black, strokewidth=1.5,  label= L"$\Delta Pr$ [mm]")

         barplot!(Ax1, 1:Nit, - Df[!,:ΔQ_sum], color=:red, strokecolor=:black, strokewidth=1.5,  label=L"$\Delta Q$ [mm]")

         lines!(Ax1, 1:Nit, Df[!,:ΔSink_sum], linewidth=2, colour=:green,  label= L"$\Delta EvapTransp$ [mm]")
      else
         barplot!(Ax1, 1:Nit,Df[!,:ΔPrThroughfall_sum_sum], color=:blue, strokecolor=:black, strokewidth=1.5,  label= L"$\Delta Pr$ [mm]")

         barplot!(Ax1, 1:Nit, - Df[!,:ΔQ_sum_sum], color=:red, strokecolor=:black, strokewidth=1.5,  label=L"$\Delta Q$ [mm]")

         lines!(Ax1, 1:Nit, Df[!,:ΔSink_sum_sum], linewidth=2, colour=:green,  label= L"$\Delta EvapTransp$ [mm]")
      end

         Leg = Legend(Fig[1,2], Ax1, framevisible=true, orientation=:horizontal, tellheight=true, tellwidth=true, nbanks=3, labelsize=14)
     
   Yticks = Header_θZₐᵥₑᵣ
   Ax2 = Axis(Fig[3, 1], xticks = (1:Nit, Xticks), yticks = (1:N_θZₐᵥₑᵣ, Yticks), height=Height, width=Width, ylabel= L"$Zcell$ [mm]")
      Ax2.xticklabelrotation = π / 3

      θZₐᵥₑᵣ_2D = Matrix{Float64}(Df[!, Df_θZₐᵥₑᵣ])

      θZₐᵥₑᵣ_2D = reverse(θZₐᵥₑᵣ_2D, dims=2)

      # :deepsea
      Hmap = heatmap!(Ax2, θZₐᵥₑᵣ_2D, colorrange=(0.0, 1.0), colormap = Reverse(:linear_kbc_5_95_c73_n256))

      Colorbar(Fig[3, 2], Hmap; label="Se[-]", width=20, ticks = 0:0.25:1)

   for iT in 1:Nit, iZ in 1:N_θZₐᵥₑᵣ
      txtcolor = θZₐᵥₑᵣ_2D[iT, iZ] > 0.5 ? :white : :black
      text!(Ax2, "$(round(θZₐᵥₑᵣ_2D[iT,iZ], digits = 2))", position = (iT, iZ),
         color = txtcolor, align = (:center, :center), textsize=18)
   end

   Ax2.xticklabelalign = (:right, :center)

   colgap!(Fig.layout, 10)
   rowgap!(Fig.layout, 10)
   resize_to_layout!(Fig)
   trim!(Fig.layout)
   save(PathOutput, Fig) # size = 600 x 450 pt
   display(Fig)

return nothing
end

PLOT_HEATMAP()


