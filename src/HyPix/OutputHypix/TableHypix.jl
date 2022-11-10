# =============================================================
#		module: tableHyix
# =============================================================
module tableHypix

   import ..cst, ..tool, ..wrc, ..kunsat, ..θaver
   using Dates
   import CSV, Tables
   using Statistics
   using DataFrames, CSV, Query, Statistics, Dates
   export TABLE_HYPIX

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : TABLE_HYPIX
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function TABLE_HYPIX(∑∑ΔSink, ∑Pet_Net, ∑Pr_Clim, ∑Q_Z, ∑T, ∑WaterBalanceη_Reduced, ∑ΔQ_Bot, CccBest, Date_Reduced, discret, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, hydro, hydroHorizon, i∑T_CalibrStart, iMultistep, iNonConverge_iOpt, iOpt_Count, iScenario, N_SoilLayer, Nit, Nit_Reduced, NseBest, Nz, optionHypix, paramHypix, pathInputHypix, pathOutputHypix, SiteName, veg, WilmotBest, WofBest, Z, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPrGross_Reduced, ΔPrThroughfall_Reduced, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔRunTimeHypix, ΔSink_Reduced, ΔT_Average, θ_Reduced, θroot_Mean, Ψ_Reduced; θobs_Reduced=[0], θsim_Aver=[0])

			println("		=== === START: Table === ===")

         if optionHypix.opt.Optimisation
            iSim = iMultistep
         else
            iSim = iOpt_Count
         end

         # Writing values of hydraulic parameters
         tableHypix.HYDRO(hydroHorizon, iSim, pathOutputHypix)

         # Writing values of veg parameters
         tableHypix.VEG(veg, iSim, pathOutputHypix)

         tableHypix.PERFORMANCE(∑∑ΔSink, ∑Pet_Net, ∑Pr_Clim, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, i∑T_CalibrStart, iNonConverge_iOpt, Nit, NseBest, paramHypix, pathOutputHypix, SiteName, WilmotBest, WofBest, ΔRunTimeHypix, ∑T, ΔT_Average, θroot_Mean)

            if optionHypix.special.Zreduced > 1.0
               tableHypix. PERFORMANCE_Q600(∑∑ΔSink, ∑Pet_Net, ∑Pr_Clim, ∑Q_Z, ∑T, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, i∑T_CalibrStart, iNonConverge_iOpt, Nit, NseBest, paramHypix, pathOutputHypix, SiteName, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average, θroot_Mean)	
            end

         if optionHypix.Table_Discretization
            tableHypix.DISCRETISATION_RRE(discret, Nz, Z[1:Nz], pathOutputHypix)
         end
         if optionHypix.Table_TimeSerie
            tableHypix.TIME_SERIES_DAILY( ∑WaterBalanceη_Reduced[1:Nit_Reduced], Date_Reduced[1:Nit_Reduced], iSim, iScenario, Nit_Reduced, pathInputHypix, pathOutputHypix, ΔEvaporation_Reduced[1:Nit_Reduced], ΔPet_Reduced[1:Nit_Reduced], ΔPond_Reduced[1:Nit_Reduced], ΔQ_Reduced[1:Nit_Reduced,1], ΔPrGross_Reduced[1:Nit_Reduced], ΔQ_Obs_Reduced, ΔQ_Reduced[1:Nit_Reduced,Nz+1], ΔRunoff_Reduced[1:Nit_Reduced], ΔSink_Reduced[1:Nit_Reduced])
         end
         if optionHypix.Table_θ
            tableHypix.θ(Date_Reduced[1:Nit_Reduced], θ_Reduced[1:Nit_Reduced,1:Nz], discret.Znode[1:Nz], iSim, pathOutputHypix)
         end
         if optionHypix.Table_Se
            tableHypix.SE(Date_Reduced[1:Nit_Reduced], hydro, iSim,  Nz, pathOutputHypix, discret.Znode[1:Nz], θ_Reduced[1:Nit_Reduced,1:Nz])
         end
         if optionHypix.Table_Ψ
            tableHypix.Ψ(Date_Reduced[1:Nit_Reduced], Ψ_Reduced[1:Nit_Reduced,1:Nz], discret.Znode[1:Nz], iSim, pathOutputHypix)
         end
         if optionHypix.Table_Q
            tableHypix.Q(Date_Reduced[1:Nit_Reduced], ΔQ_Reduced[1:Nit_Reduced,1:Nz+1], Z[Nz], discret.Znode[1:Nz], iSim, pathOutputHypix)
         end
         if optionHypix.Table_θΨ
            tableHypix.θΨ(hydroHorizon, iSim, N_SoilLayer, optionHypix, paramHypix, pathOutputHypix)
            tableHypix.KΨ(hydroHorizon, iSim, N_SoilLayer, optionHypix, paramHypix, pathOutputHypix)
         end
         if optionHypix.θavr_RootZone && optionHypix.θobs
            tableHypix.θAVERAGE(Date_Reduced[1:Nit_Reduced], iSim, θobs_Reduced[1:Nit_Reduced], θsim_Aver[1:Nit_Reduced], pathOutputHypix)
         end
         if optionHypix.Table_θZₐᵥₑᵣ
            tableHypix.θZAVER(Date_Reduced, discret, iSim, Nz, paramHypix, pathOutputHypix, Z, θ_Reduced)
         end
         if optionHypix.Table_Statistic ||  optionHypix.Plot_HeatMap_MonthCombine || optionHypix.Plot_HeatMap_YearMonth
            tableHypix.θSTATISTICS(Date_Reduced, discret, hydro, iSim, Nz, paramHypix, pathOutputHypix, Z, θ_Reduced, ΔPrThroughfall_Reduced[1:Nit_Reduced], ΔRunoff_Reduced, ΔQ_Reduced[1:Nit_Reduced, Nz+1], ΔSink_Reduced[1:Nit_Reduced])
         end

      println("		=== === END: Table === === \n")
      return nothing
      end  # function: TABLE_HYPIX
   # ------------------------------------------------------------------

   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : PERFORMACE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function PERFORMANCE(∑∑ΔSink, ∑Pet_Net, ∑Pr_Clim, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, i∑T_CalibrStart, iNonConverge_iOpt, Nit, NseBest, paramHypix, pathOutputHypix, SiteName, WilmotBest, WofBest, ΔRunTimeHypix, ∑T, ΔT_Average, θroot_Mean; Option_Yearly=true)	

         if Option_Yearly
            Convert_Yearly = 365.0 * cst.Day_2_Second / (∑T[Nit] - ∑T[i∑T_CalibrStart]) 
         else
            Convert_Yearly = 1.0
         end

         Date=now()
         CurrentDate ="_Year" * string(year(Date)) * "_Month" *  string(month(Date)) * "_Day" * string(day(Date)) 

         Path = pathOutputHypix.Table_Performance * CurrentDate * ".csv"
         Path = pathOutputHypix.Table_Performance * "Testing" * ".csv"

         Header = [ "Multisteps", "Wof[mm]", "NseBest[-]" ,"Ccc[-]", "Wilmot[-]", "∑Pr[mm year⁻¹]", "∑Pet[mm year⁻¹]","∑ΔSink[mm year⁻¹]","∑ΔQ_Bot[mm year⁻¹]" , "θroot_Mean[mm]" ,"Global_WaterBalance[mm]" ,"Global_WaterBalance_NormPr[-]", "ΔT_Average[second]","iNonConverge[count]" , "Efficiency[count/Hour]", "ΔRunTimeHypix[second]"] # "SiteName", 

         Id = 1:1:length(WofBest)



         CSV.write(Path, Tables.table([Id WofBest NseBest CccBest WilmotBest Convert_Yearly*∑Pr_Clim Convert_Yearly*∑Pet_Net Convert_Yearly*∑∑ΔSink Convert_Yearly*∑ΔQ_Bot θroot_Mean Global_WaterBalance Global_WaterBalance_NormPr ΔT_Average iNonConverge_iOpt Efficiency ΔRunTimeHypix]), writeheader=true, header=Header, bom=true)
      return nothing
      end # function PERFORMACE
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : PERFORMACE QZ600
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function PERFORMANCE_Q600(∑∑ΔSink, ∑Pet_Net, ∑Pr_Clim, ∑Q_Z, ∑T, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, i∑T_CalibrStart, iNonConverge_iOpt, Nit, NseBest, paramHypix, pathOutputHypix, SiteName, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average, θroot_Mean; Option_Yearly=true)

         if Option_Yearly
            Convert_Yearly = 365.0 * cst.Day_2_Second / (∑T[Nit] - ∑T[i∑T_CalibrStart]) 
         else
            Convert_Yearly = 1.0
         end
         
         Path = pathOutputHypix.Table_Performance * "_Q600" * ".csv"

         Header = ["SiteName" , "Multisteps", "Wof[mm]", "NseBest[-]" ,"Ccc[-]", "Wilmot[-]", "∑Pr[mm year⁻¹]", "∑Pet[mm year⁻¹]" ,"∑ΔSink[mm year⁻¹]","∑ΔQ_Bot[mm year⁻¹]", "∑ΔQ600[mm year⁻¹]" , "Error600[%]","θroot_Mean[mm]" ,"Global_WaterBalance[mm]" ,"Global_WaterBalance_NormPr[-]", "ΔT_Average[second]","iNonConverge[count]" , "Efficiency[count/Hour]", "ΔRunTimeHypix[second]"]
      
         Id = 1:1:length(WofBest)

         CSV.write(Path, Tables.table([SiteName Id WofBest NseBest CccBest WilmotBest Convert_Yearly*∑Pr_Clim Convert_Yearly*∑Pet_Net Convert_Yearly*∑∑ΔSink Convert_Yearly*∑ΔQ_Bot Convert_Yearly*∑Q_Z 100.0.*abs.(∑ΔQ_Bot.-∑Q_Z)./∑ΔQ_Bot θroot_Mean Global_WaterBalance Global_WaterBalance_NormPr ΔT_Average iNonConverge_iOpt Efficiency ΔRunTimeHypix]), writeheader=true, header=Header, bom=true)
       return nothing
      end # function PERFORMACE
   #------------------------------------------------------
   

   # ===================================================
   #          TimeStep daily
   # ===================================================
      function TIME_SERIES_DAILY(∑WaterBalanceη_Reduced, Date_Reduced, iSim, iScenario, Nit_Reduced, pathInputHypix, pathOutputHypix, ΔEvaporation_Reduced, ΔPet_Reduced, ΔPond_Reduced, ΔPr_Soil, ΔPrGross_Reduce, ΔQ_Obs_Reduced, ΔQ_Reduced, ΔRunoff_Reduced, ΔSink_Reduced)

         Path = pathOutputHypix.Table_TimeSerie_Daily * "_" * string(iSim) * ".csv"

         Id = 1:1:Nit_Reduced

         Year₁   = fill(0::Int64, Nit_Reduced)
         Month₁  = fill(0::Int64, Nit_Reduced)
         Day₁    = fill(0::Int64, Nit_Reduced)
         Hour₁   = fill(0::Int64, Nit_Reduced)
         Minute₁ = fill(0::Int64, Nit_Reduced)
         Second₁ = fill(0::Int64, Nit_Reduced)

         for iT=1:Nit_Reduced
            Year₁[iT]   = year(Date_Reduced[iT])
            Month₁[iT]  = month(Date_Reduced[iT])
            Day₁[iT]    = day(Date_Reduced[iT])
            Hour₁[iT]   = hour(Date_Reduced[iT])
            Minute₁[iT] = minute(Date_Reduced[iT])
            Second₁[iT] = second(Date_Reduced[iT])
         end

         if !(isempty(pathInputHypix.Drainage[iScenario]))

            Header =  ["iD" ,"Year" ,"Month" ,"Day" ,"Hour" ,"Minute" ,"Second" ,"ΔPrGross[mm]", "ΔPrSoil[mm]" , "ΔPet[mm]" ,"ΔSink[mm]","ΔTranspiration[mm]" ,"ΔEvaporation[mm]" ,"ΔDrainage[mm]" ,"Hpond[mm]" ,"ΔRunoff[mm]","∑WaterBalance_η[mm]", "ΔDrainage_Obs[mm]"]

            CSV.write(Path, Tables.table([Id Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ ΔPrGross_Reduce ΔPr_Soil ΔPet_Reduced ΔSink_Reduced ΔSink_Reduced.-ΔEvaporation_Reduced ΔEvaporation_Reduced ΔQ_Reduced ΔPond_Reduced ΔRunoff_Reduced ∑WaterBalanceη_Reduced ΔQ_Obs_Reduced]), writeheader=true, header=Header, bom=true)
         else
            Header =  ["iD" ,"Year" ,"Month" ,"Day" ,"Hour" ,"Minute" ,"Second" ,"ΔPrGross[mm]", "ΔPrSoil[mm]" , "ΔPet[mm]" ,"ΔSink[mm]","ΔTranspiration[mm]" ,"ΔEvaporation[mm]" ,"ΔRecharge[mm]" ,"Hpond[mm]" ,"ΔRunoff[mm]","∑WaterBalance_η[mm]"]

            CSV.write(Path, Tables.table([Id Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ ΔPrGross_Reduce ΔPr_Soil ΔPet_Reduced ΔSink_Reduced ΔSink_Reduced.-ΔEvaporation_Reduced ΔEvaporation_Reduced ΔQ_Reduced ΔPond_Reduced ΔRunoff_Reduced ∑WaterBalanceη_Reduced]), writeheader=true, header=Header, bom=true)
         end
   
      return nothing
      end # Table  TIME_SERIES_DAILY
   #------------------------------------------------------


   # ===================================================
   #          θ
   # ===================================================
      function θ(Date_Reduced, θ_Reduced, Znode, iSim, pathOutputHypix)
         Path = pathOutputHypix.Table_θ * "_" * string(iSim) * ".csv"

         Nit_Reduced = length(Date_Reduced)

         Year₁   = fill(0::Int64, Nit_Reduced)
         Month₁  = fill(0::Int64, Nit_Reduced)
         Day₁    = fill(0::Int64, Nit_Reduced)
         Hour₁   = fill(0::Int64, Nit_Reduced)
         Minute₁ = fill(0::Int64, Nit_Reduced)
         Second₁ = fill(0::Int64, Nit_Reduced)

         for iT=1:Nit_Reduced
            Year₁[iT]   = year(Date_Reduced[iT])
            Month₁[iT]  = month(Date_Reduced[iT])
            Day₁[iT]    = day(Date_Reduced[iT])
            Hour₁[iT]   = hour(Date_Reduced[iT])
            Minute₁[iT] = minute(Date_Reduced[iT])
            Second₁[iT] = second(Date_Reduced[iT])
         end

         # Adding an other column
            Header = ["θ[mm² mm⁻²]_Year", "Month", "Day", "Hour", "Minute", "Second Znode[mm]"]
            Header = vcat(Header, string.(-Znode))
      
         CSV.write(Path, Tables.table([Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ θ_Reduced]), writeheader=true, header=Header, bom=true)
      return nothing
      end  # Table θ
   #------------------------------------------------------


   # ===================================================
   #          SE
   # ===================================================
      function SE(Date_Reduced, hydro, iSim, Nz, pathOutputHypix, Znode, θ_Reduced; 🏁Θr=true)
       
         Path = pathOutputHypix.Table_Se * "_" * string(iSim) * ".csv"

         Nit_Reduced = length(Date_Reduced)

         Year₁   = fill(0::Int64, Nit_Reduced)
         Month₁  = fill(0::Int64, Nit_Reduced)
         Day₁    = fill(0::Int64, Nit_Reduced)
         Hour₁   = fill(0::Int64, Nit_Reduced)
         Minute₁ = fill(0::Int64, Nit_Reduced)
         Second₁ = fill(0::Int64, Nit_Reduced)

         Se₂ = fill(0.0::Float64, Nit_Reduced, Nz)

         for iT=1:Nit_Reduced
            Year₁[iT]   = (year(Date_Reduced[iT]))
            Month₁[iT]  = (month(Date_Reduced[iT]))
            Day₁[iT]    = (day(Date_Reduced[iT]))
            Hour₁[iT]   = (hour(Date_Reduced[iT]))
            Minute₁[iT] = (minute(Date_Reduced[iT]))
            Second₁[iT] = (second(Date_Reduced[iT]))

            for iT=1:Nit_Reduced, iZ=1:Nz
               if  🏁Θr
                  Se₂[iT, iZ] = (θ_Reduced[iT, iZ] - hydro.θr[iZ]) / (hydro.θs[iZ] - hydro.θr[iZ])
               else
                  Se₂[iT, iZ] = θ_Reduced[iT, iZ] / hydro.θs[iZ]
               end
            end
         end
         
          # Adding an other column
            Header = ["Se[0-1]_Year", "Month", "Day", "Hour", "Minute", "Second Znode[mm]"]
            Header = vcat(Header, string.(-Znode))
    
          CSV.write(Path, Tables.table([Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ Se₂]), writeheader=true, header=Header, bom=true)
      return nothing
      end  # Table θ
   #------------------------------------------------------


   # ===================================================
   #          θZAVER
   # ===================================================
      function θZAVER(Date_Reduced, discret, iSim, Nz, paramHypix, pathOutputHypix, Z, θ_Reduced)
         
         Path = pathOutputHypix.Table_θZₐᵥₑᵣ * "_" * string(iSim) * ".csv"

         Nit_Reduced = length(Date_Reduced)

         Year₁   = fill(0::Int64, Nit_Reduced)
         Month₁  = fill(0::Int64, Nit_Reduced)
         Day₁    = fill(0::Int64, Nit_Reduced)
         Hour₁   = fill(0::Int64, Nit_Reduced)
         Minute₁ = fill(0::Int64, Nit_Reduced)
         Second₁ = fill(0::Int64, Nit_Reduced)

         for iT=1:Nit_Reduced
            Year₁[iT]   = year(Date_Reduced[iT])
            Month₁[iT]  = month(Date_Reduced[iT])
            Day₁[iT]    = day(Date_Reduced[iT])
            Hour₁[iT]   = hour(Date_Reduced[iT])
            Minute₁[iT] = minute(Date_Reduced[iT])
            Second₁[iT] = second(Date_Reduced[iT])
         end

         # Computing θ average at the incremented depths described in paramHypix.table.θZₐᵥₑᵣ
            θZₐᵥₑᵣ = θaver.θAVER(discret; Z=Z, θ_Reduced=θ_Reduced, Nz=Nz, Nit_Reduced=Nit_Reduced, ZaverArray=paramHypix.table.θZₐᵥₑᵣ)
         
         Header = ["Date", "Year", "Month", "Day", "Hour", "Minute", "Second"]
         Header_θZₐᵥₑᵣ =  "Z" .* string.(paramHypix.table.θZₐᵥₑᵣ) .*"mm"
         Header = vcat(Header,  Header_θZₐᵥₑᵣ)
      
         CSV.write(Path, Tables.table([Date_Reduced year.(Date_Reduced[:])  month.(Date_Reduced[:]) day.(Date_Reduced[:]) hour.(Date_Reduced[:])  minute.(Date_Reduced[:]) second.(Date_Reduced[:]) θZₐᵥₑᵣ]), writeheader=true, header=Header, bom=true)
      return nothing
      end  # Table θ
   #------------------------------------------------------


   # ===================================================
   #          θSTATISTICS
   # ===================================================
      function θSTATISTICS(Date_Reduced, discret, hydro, iSim, Nz, paramHypix, pathOutputHypix, Z, θ_Reduced, ΔPrThroughfall, ΔRunoff, ΔQ, ΔSink; 🏁Θr=false, Nmonths_Combine=3)

         Nit_Reduced = length(Date_Reduced)
       
         # AVERAGE OVER YEAR MONTH FOR θ =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         # Computing θ average at the incremented depths described in paramHypix.table.θZₐᵥₑᵣ
            θZₐᵥₑᵣ = θaver.θAVER(discret; Z=Z, θ_Reduced=θ_Reduced, Nz=Nz, Nit_Reduced=Nit_Reduced, ZaverArray=paramHypix.table.θZₐᵥₑᵣ)
            
            Header_θZₐᵥₑᵣ =  "Z" .* string.(paramHypix.table.θZₐᵥₑᵣ) .*"mm"

            Df  = DataFrame(θZₐᵥₑᵣ, Header_θZₐᵥₑᵣ)
            Df2 = DataFrame(Date=Date_Reduced)
            Df  = hcat(Df, Df2)
            Df3 = DataFrame(ΔPrThroughfall=ΔPrThroughfall, ΔQ=ΔQ, ΔSink=ΔSink, ΔRunoff=ΔRunoff)
            Df  = hcat(Df, Df3)

            Df_YearMonth = combine(groupby(transform(Df, :Date => ByRow(yearmonth)), :Date_yearmonth), Header_θZₐᵥₑᵣ .=> mean, :ΔPrThroughfall .=> sum, :ΔQ .=> sum, :ΔSink .=> sum, :ΔRunoff .=> sum)

            Df_YearMonth[!,:Year] = map((X)->X[1] , Df_YearMonth[!,:Date_yearmonth])
            Df_YearMonth[!,:Month] = map((X)->X[2] , Df_YearMonth[!,:Date_yearmonth])
            
            # Removing :Date_yearmonth
            select!(Df_YearMonth, Not(:Date_yearmonth))

            # Reordering
            select!(Df_YearMonth, :Year, :Month, Not([:Year, :Month]))

            Path = pathOutputHypix.Table_Statistic * "_θZaver_YearMonth_" * string(iSim) * ".csv"
            CSV.write(Path, Df_YearMonth, bom=true)


         # AVERAGE OVER YEAR MONTH FOR SE =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            Se₂ = fill(0.0::Float64, Nit_Reduced, Nz)
            for iT=1:Nit_Reduced
               for iT=1:Nit_Reduced, iZ=1:Nz
                  if  🏁Θr
                     Se₂[iT, iZ] = (θ_Reduced[iT, iZ] - hydro.θr[iZ]) / (hydro.θs[iZ] - hydro.θr[iZ])
                  else
                     Se₂[iT, iZ] = θ_Reduced[iT, iZ] / hydro.θs[iZ]
                  end
               end
            end
            SeZₐᵥₑᵣ = θaver.θAVER(discret; Z=Z, θ_Reduced=Se₂, Nz=Nz, Nit_Reduced=Nit_Reduced, ZaverArray=paramHypix.table.θZₐᵥₑᵣ)

            Df  = DataFrame(SeZₐᵥₑᵣ, Header_θZₐᵥₑᵣ)
            Df2 = DataFrame(Date=Date_Reduced)
            Df  = hcat(Df, Df2)
            Df3 = DataFrame(ΔPrThroughfall=ΔPrThroughfall, ΔQ=ΔQ, ΔSink=ΔSink, ΔRunoff=ΔRunoff )
            Df  = hcat(Df, Df3)

            Df_YearMonth = combine(groupby(transform(Df, :Date => ByRow(yearmonth)), :Date_yearmonth), Header_θZₐᵥₑᵣ .=> mean, :ΔPrThroughfall .=> sum, :ΔQ .=> sum, :ΔSink .=> sum, :ΔRunoff .=> sum)
            Df_YearMonth[!,:Year] = map((X)->X[1] , Df_YearMonth[!,:Date_yearmonth])
            Df_YearMonth[!,:Month] = map((X)->X[2] , Df_YearMonth[!,:Date_yearmonth])
            
            # Removing :Date_yearmonth
            select!(Df_YearMonth, Not(:Date_yearmonth))

            # Reordering
            select!(Df_YearMonth, :Year, :Month, Not([:Year, :Month]))

            Path = pathOutputHypix.Table_Statistic * "_SeZaver_YearMonth_" * string(iSim) * ".csv"
            CSV.write(Path, Df_YearMonth, bom=true)

         
         # AVERAGE OVER YEAR AND 3 MONTH FOR SE =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         # Updated names of Header
            Ncombine = length(Df_YearMonth[!,:Month])
            Header_θZₐᵥₑᵣ = Header_θZₐᵥₑᵣ .* "_mean"
            MonthlyCombine = fill(""::String, Ncombine)

            for iT=1:Ncombine
               iMonthsCombine = max(floor(Df_YearMonth[iT,:Month] / 3), 1)
               iMonth_End = Int64(iMonthsCombine * Nmonths_Combine)
               iMonth_Start = Int64(iMonth_End - (Nmonths_Combine-1))

               MonthlyCombine[iT] = string(Df_YearMonth[iT,:Year]) * "_Month_" * string(iMonth_Start) * "_"  * string(iMonth_End)
            end

            Df = DataFrame(MonthlyCombine=MonthlyCombine)

            Df_YearMonth  = hcat(Df_YearMonth, Df)

            Df_YearMonthCombine = combine(groupby(Df_YearMonth, :MonthlyCombine), Header_θZₐᵥₑᵣ .=> mean, :ΔPrThroughfall_sum .=> sum, :ΔQ_sum .=> sum, :ΔSink_sum .=> sum, :ΔRunoff_sum .=> sum)

            rename!(Df_YearMonthCombine,:ΔPrThroughfall_sum_sum => :ΔPrThroughfall_sum)
            rename!(Df_YearMonthCombine,:ΔQ_sum_sum => :ΔQ_sum)
            rename!(Df_YearMonthCombine,:ΔSink_sum_sum => :ΔSink_sum)
            rename!(Df_YearMonthCombine,:ΔRunoff_sum_sum => :ΔRunoff_sum)
            Header_θZₐᵥₑᵣ_Mean = Header_θZₐᵥₑᵣ .* "_mean"
            rename!(Df_YearMonthCombine, Header_θZₐᵥₑᵣ_Mean .=> Header_θZₐᵥₑᵣ)

            N =length(Df_YearMonthCombine[!, :MonthlyCombine])

            Year₁ = fill("", N)
            Month₁ = fill("", N)
            MonthlyCombine_Data = Vector{String}(Df_YearMonthCombine[!, :MonthlyCombine])

            for (i, iMonthlyCombine) in enumerate(MonthlyCombine_Data)
               Year₁[i] , ~ , Month₁[i], ~ = split(iMonthlyCombine,"_")
            end

            Year₁ = parse.(Int64, Year₁)
            Month₁ =parse.(Int64, Month₁)

            Df_YearMonthCombine[!,:Year] = Year₁
            Df_YearMonthCombine[!,:Month] = Month₁

            Path = pathOutputHypix.Table_Statistic * "_SeZaver_YearMonth_Combine_" * string(iSim) * ".csv"
            CSV.write(Path, Df_YearMonthCombine, bom=true)

          
         # AVERAGE OVER WATER YEAR FOR FLUXES =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            # WATER YEAR
            # Water year in New Zealand is from the 1rst of July and ends on the 30th of June: 
            Date_WaterYear = fill(Date_Reduced[1]::Dates.DateTime, Nit_Reduced)
               for iT=1:Nit_Reduced
                  Date_WaterYear[iT] = (Date_Reduced[iT] - Dates.Day(181))
               end
               Df = DataFrame(Date=Date_WaterYear, ΔPrThroughfall=ΔPrThroughfall, ΔQ=ΔQ, ΔSink=ΔSink )

               Df_Year = combine(groupby(transform(Df, :Date => ByRow(year)), :Date_year), :ΔPrThroughfall .=> sum, :ΔQ .=> sum, :ΔSink .=> sum)

               rename!(Df_Year,:Date_year => :WaterYear)

               Path = pathOutputHypix.Table_Statistic * "_SeZaver_Year_" * string(iSim) * ".csv"
               CSV.write(Path, Df_Year, bom=true)
      return nothing
      end  # Table θ
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θAVERAGE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function θAVERAGE(Date_Reduced, iSim, θobs_Reduced, θsim_Aver, pathOutputHypix)

         Path = pathOutputHypix.Table_θaverage * string(iSim) * ".csv"

         Header = ["Id", "Year","Month","Day" ,"ThetaObs_Aver", "ThetaSim_Aver"]

         Id = 1:1:length(θsim_Aver)
         Year = year.(Date_Reduced)
         Month = month.(Date_Reduced)
         Day = day.(Date_Reduced)

         CSV.write(Path, Tables.table([Id Year Month Day θobs_Reduced θsim_Aver]), writeheader=true, header=Header, bom=true)
      return nothing			
      end # function: θAVERAGE
   #------------------------------------------------------


   # ===================================================
   #          Q
   # ===================================================
      function Q(Date_Reduced, ΔQ_Reduced, Z_Bottom, Znode, iSim, pathOutputHypix)

         Path = pathOutputHypix.Table_Q * "_" * string(iSim) * ".csv"

         Nit_Reduced = length(Date_Reduced)

         Year₁   = fill(0::Int64, Nit_Reduced)
         Month₁  = fill(0::Int64, Nit_Reduced)
         Day₁    = fill(0::Int64, Nit_Reduced)
         Hour₁   = fill(0::Int64, Nit_Reduced)
         Minute₁ = fill(0::Int64, Nit_Reduced)
         Second₁ = fill(0::Int64, Nit_Reduced)

         for iT=1:Nit_Reduced
            Year₁[iT]   = year(Date_Reduced[iT])
            Month₁[iT]  = month(Date_Reduced[iT])
            Day₁[iT]    = day(Date_Reduced[iT])
            Hour₁[iT]   = hour(Date_Reduced[iT])
            Minute₁[iT] = minute(Date_Reduced[iT])
            Second₁[iT] = second(Date_Reduced[iT])
         end
         
         # Adding an other column
         append!(Znode, Z_Bottom)

         Header = ["ΔFlux[mm]_Year", "Month", "Day", "Hour", "Minute", "Second Znode[mm]"]
         Header = vcat(Header, string.(-Znode))

         CSV.write(Path, Tables.table([Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ ΔQ_Reduced]), writeheader=true, header=Header, bom=true) 
      return nothing
      end  # function Q
   #------------------------------------------------------


   # ===================================================
   #          Ψ
   # ===================================================
      function Ψ(Date_Reduced, Ψ_Reduced, Znode, iSim, pathOutputHypix)

         Path = pathOutputHypix.Table_Ψ * "_" * string(iSim) * ".csv"

         Nit_Reduced = length(Date_Reduced)

         Year₁   = fill(0::Int64, Nit_Reduced)
         Month₁  = fill(0::Int64, Nit_Reduced)
         Day₁    = fill(0::Int64, Nit_Reduced)
         Hour₁   = fill(0::Int64, Nit_Reduced)
         Minute₁ = fill(0::Int64, Nit_Reduced)
         Second₁ = fill(0::Int64, Nit_Reduced)

         for iT=1:Nit_Reduced
            Year₁[iT]   = year(Date_Reduced[iT])
            Month₁[iT]  = month(Date_Reduced[iT])
            Day₁[iT]    = day(Date_Reduced[iT])
            Hour₁[iT]   = hour(Date_Reduced[iT])
            Minute₁[iT] = minute(Date_Reduced[iT])
            Second₁[iT] = second(Date_Reduced[iT])
         end

         Header = ["Ψ[mm] Year", "Month", "Day", "Hour", "Minute", "Second Znode[mm"]
         Header = vcat(Header, string.(-Znode))

         CSV.write(Path, Tables.table([Year₁ Month₁ Day₁ Hour₁ Minute₁ Second₁ -Ψ_Reduced]), writeheader=true, header=Header, bom=true)
      return nothing
      end  # function Ψ
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θΨ
   # 		Tabular values of the hydroParam model
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function θΨ(hydroHorizon, iSim, N_SoilLayer, optionₘ, paramHypix, pathOutputHypix)

         Path = pathOutputHypix.Table_θΨ * "_" * string(iSim) * ".csv"

         N_θΨobs = Int64(length(paramHypix.ploting.θΨ_Table))

         # Writting the Header
            FieldName_String = fill(""::String, N_θΨobs)

            for i =1:N_θΨobs
               FieldName_String[i] = string(-paramHypix.ploting.θΨ_Table[i] )
            end
            pushfirst!(FieldName_String, string("Ψ[mm] / θ(Ψ)[mm² mm⁻²]")) # Write the "Id" at the very begenning
         
         # Computing θ at required θ
            θ_Mod = fill(0.0::Float64, (N_SoilLayer, N_θΨobs))
            for iZ=1:N_SoilLayer, iΨ =1:N_θΨobs
                  Ψ_Mod =paramHypix.ploting.θΨ_Table[iΨ]
                  θ_Mod[iZ, iΨ] = wrc.Ψ_2_θDual(optionₘ, Ψ_Mod, iZ, hydroHorizon)
            end # iZ

         # Concatenating the 2 matrices
         Id = 1:1:N_SoilLayer

         θ_Mod = hcat(Id, θ_Mod)

         # Writting the table
            CSV.write(Path, Tables.table(θ_Mod[1:N_SoilLayer, 1:N_θΨobs+1]), writeheader=true, header=FieldName_String, bom=true)

      return nothing	
      end  # function:  θΨK_PSD
   #------------------------------------------------------

   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : KΨ
   # 		Tabular values of the hydroParam model
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function KΨ(hydroHorizon, iSim, N_SoilLayer, optionₘ, paramHypix, pathOutputHypix)

         Path = pathOutputHypix.Table_KΨ * "_" * string(iSim) * ".csv"

         N_θΨobs = Int64(length(paramHypix.ploting.θΨ_Table))

         # Writting the Header
            FieldName_String = fill(""::String, N_θΨobs)

            for i =1:N_θΨobs
               FieldName_String[i] = string(paramHypix.ploting.θΨ_Table[i])
            end
            pushfirst!(FieldName_String, string("Ψ[mm] / K(Ψ)[mm/hour]")) # Write the "Id" at the very begenning
         
         # Computing θ at required θ
            K_Mod = fill(0.0::Float64, (N_SoilLayer, N_θΨobs))
            for iZ=1:N_SoilLayer, iΨ =1:N_θΨobs
               Ψ_Mod =paramHypix.ploting.θΨ_Table[iΨ]
               K_Mod[iZ, iΨ] = kunsat.Ψ_2_KUNSAT(optionₘ, Ψ_Mod, iZ, hydroHorizon) .* cst.MmS_2_MmH
            end # iZ

         # Concatenating the 2 matrices
         Id = 1:1:N_SoilLayer

         K_Mod = hcat(Id, K_Mod)

         # Writting the table
            CSV.write(Path, Tables.table(K_Mod[1:N_SoilLayer,1:N_θΨobs+1]), writeheader=true, header=FieldName_String, bom=true)
      return nothing	
      end  # function:  θΨK_PSD
   #------------------------------------------------------


   # ===================================================
   #          DISCRETISATION AUTO
   # ===================================================
      function DISCRETISATION_AUTO(Flag_θΨini::Symbol, Layer::Vector{Float64}, PathDiscretisation::String, Z::Vector{Float64}, θini_or_Ψini_Cell::Vector{Float64})

         if Flag_θΨini == :Ψini
            Header = ["iZ","Z", "Layer", "Ψini"]

         elseif Flag_θΨini == :θini
            Header = ["iZ","Z", "Layer", "θini"]
         end

         iZ = collect(1:1:length(Z))

         CSV.write(PathDiscretisation, Tables.table([iZ Z Layer θini_or_Ψini_Cell]), writeheader=true, header=Header, bom=true)
      return nothing
      end # Table DISCRETISATION_AUTO
   #------------------------------------------------------


   # ===================================================
   #          Discretization
   # ===================================================
      function DISCRETISATION_RRE(discret, Nz, Z, pathOutputHypix)
         Header =  ["Z", "ΔZ", "ΔZ_⬓", "Znode", "ΔZ_Aver", "ΔZ_W", "Z_CellUp"]

         CSV.write(pathOutputHypix.Table_Discretisation, Tables.table( [Z[1:Nz] discret.ΔZ[1:Nz] discret.ΔZ_⬓[1:Nz] discret.Znode[1:Nz] discret.ΔZ_Aver[1:Nz] discret.ΔZ_W[1:Nz] discret.Z_CellUp[1:Nz]]), writeheader=true, header=Header, bom=true)
      return nothing
      end # Table DISCRETISATION
   #------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : HYDRO
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function HYDRO(hydroHorizon, iSim, pathOutputHypix)

         Path = pathOutputHypix.Table_Hydro  * "_" * string(iSim) * ".csv"

         N_SoilLayer = length(hydroHorizon.θs)

         Id = 1:1:N_SoilLayer

         Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilLayer, hydroHorizon)
               
         pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

         CSV.write(Path, Tables.table([Int64.(Id) Matrix]), writeheader=true, header=FieldName_String, bom=true)
      return nothing			
      end  # function: HYDRO
   #------------------------------------------------------
   

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : veg
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function VEG(veg, iSim, pathOutputHypix)
         
         Path = pathOutputHypix.Table_Veg * "_" * string(iSim) * ".csv"

         Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(1, veg)
       
         CSV.write(Path, Tables.table(Matrix), writeheader=true, header=FieldName_String, bom=true)
      return nothing
      end  # function: VEG
   #------------------------------------------------------

end  # module: tableHyix
# ............................................................