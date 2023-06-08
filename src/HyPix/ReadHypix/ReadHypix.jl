# =============================================================
#		module: readHypix
# =============================================================
module readHypix

   import ..climate, ..discretisation, ..horizonLayer, ..hydroStruct, ..memory, ..optionsHypix, ..paramsHypix, ..pathsHypix, ..tableHypix, ..thetaObs, ..tool, ..vegStruct, ..datesHypix, ..interpolate, ..hydroSmooth, ..cst
   import Dates: value, DateTime, hour, minute, month, now, Hour
   import DelimitedFiles
   import CSV, Tables
   using Base

   export READ_START, HYPIX_PARAM_OPT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : READ_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function READ_START(dateHypix, Id, iScenario, N_Scenario, Path_Hypix, pathInputHypix, ProjectHypix, SiteName)

         # READING HYPIX PARAMETERS
            # Reading HyPix parameters
               paramHypix = paramsHypix.PARAM_HYPIX(pathInputHypix.ParamHypix[iScenario])

            # Reading HyPix options
               optionHypix = optionsHypix.OPTION_HYPIX(pathInputHypix.OptionHypix[iScenario])

            # Reading HyPix paths
               pathOutputHypix = pathsHypix.PATH_HYPIX(Path_Hypix, pathInputHypix.PathHypix[iScenario], ProjectHypix, SiteName[iScenario])

         # UPDATING DATES WITH iScenario
            dateHypix = datesHypix.DATE_HYPIX(dateHypix, iScenario)
         
         # READING CLIMATE DATA ====
            clim, dateHypix = readHypix.CLIMATE(dateHypix, iScenario, pathInputHypix.Climate[iScenario])
            # Process climate data
               ∑Pet_Climate, ∑PrThroughfall_Climate, ∑T_Climate, N_∑T_Climate, Temp = climate.CLIMATE(clim, optionHypix)

         # READING DISCRETISATION ===
            # Create Discretisation.csv from SoilLayer.csv
            if optionHypix.Discretisation_File_Auto⍰ == "Auto"

               # Read SoilLayer, could be either θini, Ψini
                  🎏_θΨini, Layer, N_SoilLayer, ~, Zlayer, θini_or_Ψini = readHypix.DISCRETISATION(pathInputHypix.SoilLayer[iScenario])

               # Option smoothening hydraulic parameters
                  if optionHypix.HydroSmooth
                     Layer, N_Layer, Zlayer, θini_or_Ψini = hydroSmooth.DISCRETISATION_SMOOTENING!(🎏_θΨini, iScenario, N_SoilLayer, optionHypix, paramHypix, pathInputHypix, Zlayer, θini_or_Ψini)
                  else
                     N_Layer = copy(N_SoilLayer)
                  end

               # Performing auto discretisation			
                  Layer, Nz, Z, θini_or_Ψini = discretisation.DISCRETISATION_AUTO(Layer, optionHypix, paramHypix; N_Layer=N_Layer, Zlayer=Zlayer, θini_or_Ψini=θini_or_Ψini)

               if optionHypix.special.Zreduced > 1.0
                  Layer, Nz, Z, θini_or_Ψini = readHypix.DISCRETISATION_REDUCED(Layer, Nz, Z, θini_or_Ψini)
               end

               tableHypix.DISCRETISATION_AUTO(🎏_θΨini, Layer, pathInputHypix.Discretisation[iScenario], Z, θini_or_Ψini)
            else
               # Read discretisation
               🎏_θΨini, Layer, N_SoilLayer, Nz, Z, θini_or_Ψini = readHypix.DISCRETISATION(pathInputHypix.Discretisation[iScenario])
               N_Layer = copy(N_SoilLayer)
            end # if optionHypix.Discretisation_File_Auto⍰ == "Auto" 

            # Process discretisation of the soil profile ~~~~~
               discret = discretisation.DISCRETISATION(Nz, Z)

         # READING TABULAR CUMULATIVE ROOTDENSITY FUNCTION IF REQUIRED ===
            if !(isempty(pathInputHypix.∑RootDensityFunc[iScenario])) &&  optionHypix.RootWaterUptake
               N_iRoot, ΔRootDensity = readHypix.∑ROOTDENSITY(Nz, pathInputHypix.∑RootDensityFunc[iScenario], Z)
            else
               ΔRootDensity = empty
               N_iRoot = 0
            end

         # READING OBSERVED θ IF REQUIRED ===
            if !(isempty(pathInputHypix.θdata[iScenario]))
               optionHypix.θobs = true
               # Read observed θ
               dateHypix, obsθ = readHypix.θDATA(clim, dateHypix, iScenario, pathInputHypix.θdata[iScenario])
               # Process observed θ
                  obsθ = thetaObs.ΘOBS(obsθ, clim, discret, Z)
            else
               obsθ = []
               optionHypix.θobs = false
            end #  optionHypix.θobs

         # COMPUTING Δ∑T
            dateHypix = datesHypix.DATE_HYPIX_Δ∑T(clim, dateHypix)
      
         # LOOKUP TABLE ===
            # Initialiozing vegetation parameters into veg structure
               veg = vegStruct.VEGSTRUCT() 

               # if optionHypix.LookupTable_Lai
               #    Laiᵀ_η = readHypix.LOOKUPTABLE_LAI(clim, optionHypix, pathInputHypix.LookUpTable_Lai[iScenario], veg)
               # else
               #    Laiᵀ_η = []
               # end

               if optionHypix.LookUpTable_CropCoeficient
                  CropCoeficientᵀ_η = readHypix.LOOKUPTABLE_CROPCOEFICIENT(clim, optionHypix, pathInputHypix.LookUpTable_Crop[iScenario], veg)
               else
                  CropCoeficientᵀ_η = Float64[]
               end
               
         # INITIALIZING THE STRUCTURE ===
            # Initializing hydroHorizon structure
               hydroHorizon = hydroStruct.HYDROSTRUCT(optionHypix, N_Layer)

            # Optimisation
               hydroHorizon_best = hydroStruct.HYDROSTRUCT(optionHypix, N_Layer)
               hydro_best        = hydroStruct.HYDROSTRUCT(optionHypix, Nz)
               veg_best          = vegStruct.VEGSTRUCT()

         # HYDRAULIC VEGETATION PARAMETERS
            if !(optionHypix.opt.Optimisation)
               veg, ~ = tool.readWrite.READ_STRUCT_SIMPLE(veg, pathInputHypix.Vegetation[iScenario])

               # Hydraulic parameters
                  hydroHorizon, ~ = tool.readWrite.READ_STRUCT_SIMPLE(hydroHorizon, pathInputHypix.HydroInput[iScenario])
            else
               hydroHorizon = []      
            end # optionHypix.Optimisation

            # if 🎏_ImpermeableLayer
            #    hydro.Ks[Nz] = 1.15741E-05 #[mm s-1]  
            # end  # if: 🎏_ImpermeableLayer


         # SOILPROFILE DATA
            # Reading roots
               if  !(isempty(pathInputHypix.Soilprofile[iScenario]))
                  BottomBoundary, Zroot_Max = readHypix.SOILPROFILE(pathInputHypix.Soilprofile[iScenario])
               else
                  Zroot_Max = Inf
               end

            # Reading drainage
               if  !(isempty(pathInputHypix.Drainage[iScenario]))
                  ∑ΔQ_Obs, ∑T_Qobs = readHypix.DRAINAGE(dateHypix, pathInputHypix.Drainage[iScenario])
               else
                  ∑ΔQ_Obs=[]; ∑T_Qobs=[]
               end

            # Reading LAI 
               if !(isempty(pathInputHypix.Lai[iScenario]))
                  clim = readHypix.LAI(∑T_Climate, clim, dateHypix, pathInputHypix.Lai[iScenario])
               else
                  for iT=1:clim.N_Climate
                     clim.Lai[iT] = veg.Lai
                  end
               end

         # CORRECTING ROOTING DEPTH DEPENDING ON SOIL INFORMATION
				# Maximum rooting depth should be smaller than the maximum depth of soil and the maximum soil depth given in the o folder
				   veg.Zroot = min(veg.Zroot, Z[Nz], Zroot_Max)

         # MEMORY 
            ∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑PrThroughfall, ∑T, CropCoeficientᵀ, Hpond, K_Aver_Vect, K_Aver₀_Vect, Pkₐᵥₑᵣ, Q, Residual, ΔEvaporation, ΔLnΨmax, ΔPet, ΔPrThroughfall, ΔRunoff, ΔSink, ΔT, θ, θSim, Ψ, Ψ_Min, Ψbest = memory.MEMORY(clim, N_∑T_Climate, Nz, obsθ, optionHypix, paramHypix)
         
   return ∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑PrThroughfall, ∑PrThroughfall_Climate, ∑T, ∑T_Climate, ∑T_Qobs, ∑T_Qobs, ∑ΔQ_Obs, ∑ΔQ_Obs, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, 🎏_θΨini, 🎏_θΨini, Hpond, hydro_best, hydroHorizon, hydroHorizon_best, K_Aver_Vect, K_Aver₀_Vect, Layer, N_∑T_Climate, N_iRoot, N_Layer, N_SoilLayer, Nz, obsθ, optionHypix, paramHypix, pathInputHypix, pathOutputHypix, Pkₐᵥₑᵣ, Q, Residual, Temp, veg, veg_best, Z, ΔEvaporation, ΔLnΨmax, ΔPet, ΔPrThroughfall, ΔRootDensity, ΔRunoff, ΔSink, ΔT, θ, θini_or_Ψini, θSim, Ψ, Ψ_Min, Ψbest
   end  # function: READ_START
   # ------------------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : DISCRETISATION
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function DISCRETISATION(Path::String)
         # Read data
            Data = CSV.File(Path, header=true)
            Header = string.(Tables.columnnames(Data))

            Z           = convert(Vector{Float64}, Tables.getcolumn(Data, :Z))
            Layer       = convert(Vector{Int64}, Tables.getcolumn(Data, :Layer))
            N_SoilLayer = Int64(maximum(Layer))

            Nz = length(Z)

         # Depending on the initial boundary condition 
            if "θini" ∈ Header
               θini_or_Ψini = convert(Vector{Float64}, Tables.getcolumn(Data, :θini))
               🎏_θΨini = :θini 

            elseif "Ψini" ∈ Header
               θini_or_Ψini = convert(Vector{Float64}, Tables.getcolumn(Data, :Ψini))
               🎏_θΨini = :Ψini

            else
               error("In $Path cannot find <θini> or <Ψini> in $Header")
            end
      return 🎏_θΨini, Layer, N_SoilLayer, Nz, Z, θini_or_Ψini
      end # function DISCRETISATION
   #-------------------------------------------------------------------------
 

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : CLIMATE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Base.@kwdef mutable struct CLIMATEDATA
         Date           :: Vector{DateTime}
         Pr             :: Vector{Float64}
         Pet            :: Vector{Float64}
         Temp           :: Vector{Float64}
         N_Climate      :: Int64
         Pr_Through     :: Vector{Float64}
         Lai            :: Vector{Float64}
      end

      function CLIMATE(dateHypix, iScenario::Int64, PathClimate::String; Option_ReadTemperature=false)

         # READ DATA
            Data = CSV.File(PathClimate, header=true)

            Header = string.(Tables.columnnames(Data))

            # Automatic determine the format of the climate data
               if "Year[]" ∈ Header
                  🎏_VirtualClimateStationNz = true
               else
                  🎏_VirtualClimateStationNz = false
               end

            if !🎏_VirtualClimateStationNz
               Year   = convert(Vector{Int64}, Tables.getcolumn(Data, :Year))
               N_Climate = length(Year)
               Month  = convert(Vector{Int64}, Tables.getcolumn(Data, :Month))
               Day    = convert(Vector{Int64}, Tables.getcolumn(Data,  :Day))
               Hour   = convert(Vector{Int64}, Tables.getcolumn(Data,  :Hour))
               Minute = convert(Vector{Int64}, Tables.getcolumn(Data, :Minute))
               Second = convert(Vector{Int64}, Tables.getcolumn(Data,  :Second))
               Pr     = convert(Vector{Float64}, Tables.getcolumn(Data,  Symbol("Rain(mm)")))
               Pet    = convert(Vector{Float64}, Tables.getcolumn(Data, Symbol("PET(mm)")))

               if Option_ReadTemperature 
                  Temp=  convert(Vector{Float64}, Tables.getcolumn(Data,  Symbol("Tmax(C)")))
               else
                  Temp = fill(24.0::Float64, N_Climate)
               end

            elseif  🎏_VirtualClimateStationNz
               Year   = convert(Vector{Int64}, Tables.getcolumn(Data, Symbol("Year[]")))
               N_Climate = length(Year)
               Month  = convert(Vector{Int64}, Tables.getcolumn(Data, Symbol("Month[]")))
               Day    = convert(Vector{Int64}, Tables.getcolumn(Data,  Symbol("Day[]")))
               Hour   = fill(9::Int64, N_Climate)
               Minute = fill(0::Int64, N_Climate)
               Second = fill(0::Int64, N_Climate)
               Pr     = convert(Vector{Float64}, Tables.getcolumn(Data,  Symbol("Rain[mm]")))
               Pet    = convert(Vector{Float64}, Tables.getcolumn(Data, Symbol("PET[mm]")))

               if Option_ReadTemperature 
                  Temp=  convert(Vector{Float64}, Tables.getcolumn(Data,  Symbol("TMax[oC]")))
               else
                  Temp = fill(24.0::Float64, N_Climate)
               end
            end

         # REDUCING THE NUMBER OF SIMULATIONS SUCH THAT IT IS WITHIN THE SELECTED RANGE
            Date_Start_Minimum = DateTime(Year[2], Month[2], Day[2], Hour[2], Minute[2], Second[2]) 

            Date_End_Maximum = DateTime(Year[end], Month[end], Day[end], Hour[end], Minute[end], Second[end])
                
            if !(Date_Start_Minimum ≤ dateHypix.Date_Start ≤ dateHypix.Date_SimStart ≤ dateHypix.Date_End ≤ Date_End_Maximum)
               @show Date_Start_Minimum dateHypix.Date_Start dateHypix.Date_SimStart dateHypix.Date_End Date_End_Maximum
               error("\n ***** HyPix CLIMATE Date error ∉ [Date_Start Date_SimStart Date_End] *** \n")
            end

         # SELECTING DATES OF INTEREST
            True = falses(N_Climate)
            Date = fill(Date_Start_Minimum::DateTime, N_Climate) 
            for iT=1:N_Climate
               Date[iT] = DateTime(Year[iT], Month[iT], Day[iT], Hour[iT], Minute[iT], Second[iT])

               if (dateHypix.Date_Start ≤ Date[iT] ≤ dateHypix.Date_End)
                  True[iT] = true
               end  # if: 
            end # iT=1:N_Climate

            # Need to include one date iT-1 at the beginning to compute ΔT
               iTrue_First = findfirst(True[:])
               True[iTrue_First-1] = true

            # Adjusting the starting date
                dateHypix.Date_Start = Date[iTrue_First-1]

            # New reduced number of simulations
               Date = Date[True[:]]
               Pr   = Pr[True[:]]
               Pet  = Pet[True[:]]
               Temp = Temp[True[:]]

            # Adjusting the dates
               dateHypix.Date_Start = Date[1]
               dateHypix.Date_End   = Date[end]

            # Update N_Climate	
               N_Climate = count(True[:]) # New number of data
         
         # To be used after interception model
            Pr_Through = fill(0.0::Float64, N_Climate)
            Lai        = fill(0.0::Float64, N_Climate)

      # STRUCTURE
      return CLIMATEDATA(Date, Pr, Pet, Temp, N_Climate, Pr_Through, Lai), dateHypix
      end # function: CLIMATE
   #---------------------------------------------------------------------

     
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : HYPIX_PARAM_OPT
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function HYPIX_PARAM_OPT(hydroHorizon, iMultistep::Int64, N_SoilLayer::Int64, optionHypix, paramHypix, Path::String, veg)
         # Read data
            Data = DelimitedFiles.readdlm(Path, ',')
         # Read header
            Header = Data[1,1:end]
         # Remove first READ_ROW_SELECT
            Data = Data[2:end,begin:end]

         # Readingt the type of data
            Type, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "TYPE")

         # Reading the names of the parameters
            Name, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "NAME")

            Name_Unique = unique(Name)
      
         # Reading the values of the parameters for the simulation of interest
            Param, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "SIM_$(iMultistep)")
         
         # Minimum value of the param
            Param_Min, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MIN")

         # Maximum value of the param
            Param_Max, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MAX")

         # Determening which param to optimize
            Opt, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "OPT_$(iMultistep)")

         # Maximum value of the param
            Opt_LogTransform, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "LogTransform")
            
         # Determine if we need to optimize
            if sum(Opt) ≥ 1
               🎏_Opt = true
            else
               🎏_Opt = false
            end

         """Determening if multistep optimisation is performed (not the first step)
         This is such that the optimal values of the previous optimisation step is kept in memory
         We need to determine what next param to optimize"""
            if 🎏_Opt && (iMultistep ≥ paramHypix.opt.iOptMultiStep_Start + 1)
               🎏_MultiStepOpt = true
            else
               🎏_MultiStepOpt = false 
            end
            
         # ====================================================

         # Does not matter if repeated in multistep optimisation
         ParamOpt              = []
         ParamOpt_HorizonEq    = []
         ParamOpt_Max          = []
         ParamOpt_Min          = []
         ParamOpt_Type         = []
         ParamOpt_LogTransform = []

         # Deriving hydroHorizon
         if  !(🎏_MultiStepOpt)
            hydroHorizon = hydroStruct.HYDROSTRUCT(optionHypix, N_SoilLayer)
         end

         for i in eachindex(Name_Unique)
            # Finding the position of each Param name in .csv
               indexName = findall(isequal(Name_Unique[i]), Name)

            # Values of param for every Name to put in hydroHorizon
               Param_Vect = Float64.(Param[indexName])

            if Type[i] == "hydro" && !(🎏_MultiStepOpt)
               # Putting soil param in hydroHorizon

               # θsMacMat value depends on θs
                  if Symbol(Name_Unique[i]) == :θsMacMat_ƞ
                     for iZ =1:length(Param_Vect)
                        hydroHorizon.θsMacMat[iZ] = min((hydroHorizon.θs[iZ] - hydroHorizon.θr[iZ]) * Param_Vect[iZ] + hydroHorizon.θr[iZ], hydroHorizon.θr[iZ])
                     end
                  end 
               
               setfield!(hydroHorizon, Symbol(Name_Unique[i]), Param_Vect)

               # Minimum and maximum value of the hydraulic parameters such as θs_Min and θs_Max
                  setfield!(hydroHorizon, Symbol(Name_Unique[i] * "_Min"), Float64.(Param_Min[indexName]))
                  setfield!(hydroHorizon, Symbol(Name_Unique[i] * "_Max"), Float64.(Param_Max[indexName]))

            elseif Type[i] == "veg" && !(🎏_MultiStepOpt)
               # Putting veg param in veg
               setfield!(veg, Symbol(Name_Unique[i]), Float64(Param[indexName][1]))
            end
            
            # Param to optimize. The complication is that there are different layers of hydraulic parameters which can be optimized.  
            if sum(Opt[indexName]) > 0
               # Type of parameters
                  append!(ParamOpt_Type, [Type[i]]) 

               # Appending name of param to optimize by removing dublicates
                  append!(ParamOpt, [Name_Unique[i]])

               # Appending the iHorizon of the param which will have = values. The horizon to be optimized must follow 
                  iHorizonOpt_Start = findfirst(x->x==1, Opt[indexName])
                  iHorizonOpt_End = findlast(x->x==1, Opt[indexName])

                  append!(ParamOpt_HorizonEq, [[iHorizonOpt_Start ; iHorizonOpt_End]])

               # Minimum and Maximum value of the parameter to be optimized. If we have layers than we use the value of the top layer
               iNameOpt = findfirst(x->x==Name_Unique[i], Name) + iHorizonOpt_Start - 1

               iStart = iNameOpt
               iEnd = iNameOpt + iHorizonOpt_End - iHorizonOpt_Start

               # We take the minimum to be the minimum of iHorizonOpt_Start and iHorizonOpt_End and the same for maximum
               append!(ParamOpt_Min, minimum(Param_Min[iStart:iEnd]))
               append!(ParamOpt_Max, maximum(Param_Max[iStart:iEnd]))

               # Appending name of param to perform logTransform if optimized by removing dublicates
               if sum(Opt_LogTransform[indexName]) > 0
                  append!(ParamOpt_LogTransform, [true])

                  ParamOpt_Min[end] = log1p(ParamOpt_Min[end])
                  ParamOpt_Max[end] = log1p(ParamOpt_Max[end])
               else
                  append!(ParamOpt_LogTransform, [false])
               end

               if Param_Min[iNameOpt] > Param_Max[iNameOpt]
                  error("HYPIX ERROR: $(Param_Min[iNameOpt]) < $(Name_Unique[i]) < $(Param_Max[iNameOpt]) !")
               end
            end
         end # for loop

         NparamOpt = length(ParamOpt)

         # CHECKING FOR UNCONSISTENCY WITH OPTIONS	
         if 🎏_Opt && optionHypix.opt.σ_2_Ψm⍰ ≠ "No" && "Ψm" ∈ ParamOpt
            iψm = findfirst(isequal("Ψm"), ParamOpt)[1]

            if optionHypix.opt.σ_2_Ψm⍰=="UniqueRelationship" && "Ψm" ∈ ParamOpt
               error( "**** HyPix Error: combination of options which are not possible (optionHypix.opt.σ_2_Ψm⍰==UniqueRelationship) && (Optimise=Ψm)!")

            elseif optionHypix.opt.σ_2_Ψm⍰=="Constrained" && !("Ψm" ∈ ParamOpt)
               error("*** HyPix Error: combination of options which are not possible (optionHypix.opt.σ_2_Ψm⍰==Constrained) && (not Optimising=Ψm)!")

            elseif optionHypix.opt.σ_2_Ψm⍰=="Constrained" && ParamOpt_LogTransform[iψm]==1
               error("*** optionHypix.opt.σ_2_Ψm⍰==Constrained CANNOT log transforme Ψm") 
            end
         end # 🎏_Opt

         # Putting all the parameters in  NamedTuple
         optim = (ParamOpt_Min=ParamOpt_Min, ParamOpt_Max=ParamOpt_Max, ParamOpt_HorizonEq=ParamOpt_HorizonEq, ParamOpt_Type=ParamOpt_Type, ParamOpt=ParamOpt, NparamOpt=NparamOpt, 🎏_Opt=🎏_Opt, ParamOpt_LogTransform=ParamOpt_LogTransform)

         if 🎏_Opt
            println("	=== === Optimizing the following parameters === ===")
            println("		NparamOpt=" , NparamOpt)
            println("		ParamOpt= " , optim.ParamOpt_Type .* optim.ParamOpt)
            println("		Min_Value= " , optim.ParamOpt_Min)
            println("		Max_Value= " , optim.ParamOpt_Max)
            println("		LogTransform = " , optim.ParamOpt_LogTransform)
            println("		Hydro_HorizonEq= " , optim.ParamOpt_HorizonEq)
            println("	=== === ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ === === \n")
         end
   return hydroHorizon, optim, veg
   end  # function: HYPIX_PARAM_OPT
   #---------------------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : LOOKUPTABLE
   #		Parameters as a function of time
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function LOOKUPTABLE_LAI(clim, optionHypix, PathLai::String, veg)	
         if optionHypix.LookupTable_Lai == true
            LookUpTable_Lai, ~   = tool.readWrite.READ_HEADER(PathLai, "Lai")
         end
         
         Laiᵀ_Norm = fill(0.0::Float64, clim.N_Climate) 
         for (i, Date) in enumerate(clim.Date)
            Month = month(Date)
            if optionHypix.LookupTable_Lai
               Laiᵀ_Norm[i] = LookUpTable_Lai[Month]
            else
               Laiᵀ_Norm[i] = veg.Lai
            end
         end
      return Laiᵀ_Norm
      end  # function: LOOKUPTABLE_LAI

      function LOOKUPTABLE_CROPCOEFICIENT(clim, optionHypix, PathLai::String, veg)
         if optionHypix.LookUpTable_CropCoeficient == true
            LookUpTable_CropCoeficient, ~   = tool.readWrite.READ_HEADER(PathLai, "CropCoeficient")
         end
         
         CropCoeficientᵀ_Norm = fill(0.0::Float64, clim.N_Climate) 
         for (i, Date) in enumerate(clim.Date)
            Month = month(Date)
            if optionHypix.LookUpTable_CropCoeficient == true
               CropCoeficientᵀ_Norm[i] = LookUpTable_CropCoeficient[Month]
            else
               CropCoeficientᵀ_Norm[i] = veg.CropCoeficient
            end
         end
      return CropCoeficientᵀ_Norm
      end  # function: LOOKUPTABLE_LAI
   # ------------------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : SOILPROFILE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function SOILPROFILE(PathSoilprofile)
         Data₀ = CSV.File(PathSoilprofile, header=true)
            BottomBoundary = convert(Vector{String}, Tables.getcolumn(Data₀, Symbol("BottomBoundary")))[1]
            Zroot_Max = convert(Vector{Int64}, Tables.getcolumn(Data₀, Symbol("RootingDepth[mm]")))[1]

      return BottomBoundary, Zroot_Max
      end  # function: SOILPROFILE
   # ------------------------------------------------------------------

   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : DISCRETISATION_REDUCED
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function DISCRETISATION_REDUCED(Layer, Nz, Z, θini_or_Ψini; Zmax=600.0)
         iZmax = 0
         for iZ = 1:Nz
            if Z[iZ] ≤ Zmax
               iZmax += 1
            else
               break
            end
         end # iZ=1:Nz
         Nz = iZmax
         resize!(Layer, Nz)
         resize!(Z, Nz)
         resize!(θini_or_Ψini, Nz)

      return Layer, Nz, Z, θini_or_Ψini
      end
   #---------------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : DRAINAGE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function DRAINAGE(dateHypix, Path)

         # READ DATA
            Data = CSV.File(Path, header=true)

            Year   = convert(Vector{Int64}, Tables.getcolumn(Data, :Year))
            Month  = convert(Vector{Int64}, Tables.getcolumn(Data, :Month))
            Day    = convert(Vector{Int64}, Tables.getcolumn(Data,  :Day))
            Hour   = convert(Vector{Int64}, Tables.getcolumn(Data,  :Hour))
            Minute = convert(Vector{Int64}, Tables.getcolumn(Data, :Minute))
            Second = convert(Vector{Int64}, Tables.getcolumn(Data,  :Second))

            ΔQ     = convert(Vector{Float64}, Tables.getcolumn(Data,  Symbol("Drainage[mm]")))
            N       = length(ΔQ)
         # MEMORY
            ∑ΔQ_Obs = fill(0.0::Float64, N)
            ∑T_Qobs = fill(0.0::Float64,N)
            ∑ΔQ_Obs[1] = 0.0::Float64


         # DATES &  ∑ΔQ_Obs
            Date = DateTime(Year[1], Month[1], Day[1], Hour[1], Minute[1], Second[1])
            ∑T_Qobs[1] = value(Date - dateHypix.Date_Start) / 1000
            
            for iT = 2:N
               ∑ΔQ_Obs[iT]  = ∑ΔQ_Obs[iT-1] + ΔQ[iT]

               Date = DateTime(Year[iT], Month[iT], Day[iT], Hour[iT], Minute[iT], Second[iT])
               ∑T_Qobs[iT] = value(Date - dateHypix.Date_Start) / 1000 - cst.Day_2_Second
            end

          Year=nothing; Month=nothing; Day=nothing; Hour=nothing; Minute=nothing; Second=nothing

      return ∑ΔQ_Obs, ∑T_Qobs
      end  # function: SEEPAGE
   # ------------------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : LAI
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function LAI(∑T_Climate, clim, dateHypix, Path)

         # READ DATA
            Data = CSV.File(Path, header=true)

            Year      = convert(Vector{Int64}, Tables.getcolumn(Data, :Year))
            Month     = convert(Vector{Int64}, Tables.getcolumn(Data, :Month))
            Day       = convert(Vector{Int64}, Tables.getcolumn(Data,  :Day))
            Hour      = convert(Vector{Int64}, Tables.getcolumn(Data,  :Hour))
            Minute    = convert(Vector{Int64}, Tables.getcolumn(Data, :Minute))
            Second    = convert(Vector{Int64}, Tables.getcolumn(Data,  :Second))

            Lai₀ = convert(Vector{Float64}, Tables.getcolumn(Data, :Lai))
            N    = length(Lai₀)
      
         # DATES &  ∑ΔQ_Obs
            ∑T_Lai   = fill(0.0::Float64, N)
            for iT = 1:N
               Lai_Date = DateTime(Year[iT], Month[iT], Day[iT], Hour[iT], Minute[iT], Second[iT])

               ∑T_Lai[iT] = value(Lai_Date - dateHypix.Date_Start) / 1000 - cst.Day_2_Second
            end


         # INTERPOLATING SO THAT LAI HAS THE SAME TIME STEP AS ∑T_Climate
            clim.Lai  = interpolate.INTERPOLATE_1D_LOOP(∑T_Lai, ∑T_Climate, length(∑T_Climate), N,  clim.Lai, Lai₀)


         Year=nothing; Month=nothing; Day=nothing; Hour=nothing; Minute=nothing; Second=nothing

         ∑T_Lai = nothing
      return clim
      end  # function: LAI
   # ------------------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : ∑ROOTDENSITY
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function ∑ROOTDENSITY(Nz, PathRootDensity, Z)
         Data₀ = CSV.File(PathRootDensity, header=true)

         ZrootDensity = convert(Vector{Float64}, Tables.getcolumn(Data₀, Symbol("Z[mm]")))
         ∑RootDensity = convert(Vector{Float64}, Tables.getcolumn(Data₀, Symbol("CumulRootDensityFunc[0-1]")))

         # Need to start with 0
         if ZrootDensity[1] > eps(1000.0)
            prepend!(ZrootDensity, 0.0)
            prepend!(∑RootDensity, 0.0)
         end

         Nrd = length(ZrootDensity)
         # Computing the depths of roots
            iZ = 1
            while iZ ≤ Nz && Z[iZ] ≤ ZrootDensity[end] 
               iZ +=1
            end
            N_iRoot = min(iZ - 1, Nz)
         
         ΔRootDensity     = fill(0.0::Float64, Nz)
         ∑RootDensity_Tab = fill(1.0::Float64, Nz)
         
         ∑RootDensity_Tab = interpolate.INTERPOLATE_1D_LOOP(ZrootDensity, Z, N_iRoot, Nrd, ∑RootDensity_Tab, ∑RootDensity)
         ∑RootDensity_Tab[N_iRoot] = 1.0

         ΔRootDensity[1] = ∑RootDensity_Tab[1]
         for iZ=2:Nz
            ΔRootDensity[iZ] = ∑RootDensity_Tab[iZ] - ∑RootDensity_Tab[iZ-1] 
         end
      return N_iRoot, ΔRootDensity
      end  # function: ∑ROOTDENSITY
   # ------------------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θOBSERVATION
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Base.@kwdef mutable struct θOBSERVATION
         Date      :: Vector{DateTime}
         Z  	    :: Vector{Float64}
         ithetaObs :: Vector{Int64}
         Nit       :: Int64 # Number of time steps
         Ndepth    :: Int64 # Numver of soil profile with observed θ
         θobs 	   
         ∑T  	    :: Vector{Float64}
      end # mutable struct

      function θDATA(clim, dateHypix, iScenario, Pathθobs)
      # Read data
         Data₀ = CSV.File(Pathθobs, header=true)

         Header = string.(Tables.columnnames(Data₀))

         Year   = convert(Vector{Int64}, Tables.getcolumn(Data₀, :Year))
         Month  = convert(Vector{Int64}, Tables.getcolumn(Data₀, :Month))
         Day    = convert(Vector{Int64}, Tables.getcolumn(Data₀, :Day))
         Hour   = convert(Vector{Int64}, Tables.getcolumn(Data₀, :Hour))
         Minute = convert(Vector{Int64}, Tables.getcolumn(Data₀, :Minute))
         Second = convert(Vector{Int64}, Tables.getcolumn(Data₀, :Second))

         Data = Tables.matrix(Data₀)

         Nit    = length(Year)

         # READING THE DEPTH OF Θ MEASUREMENTS FROM HEADER: data having Z=
               Array_iHeader = Int64[]
               Ndepth = 0::Int64
               iCount = 0::Int64
               for iHeader in Header
                  iCount += 1
                  if occursin("Z=", iHeader) # Searching for 'Z=' in the header
                     Ndepth += 1
                     append!(Array_iHeader, iCount) 
                  end # occursin
               end # iHeader

            # Isolating data with Z= measurements
               # Nit,~ = size(Data)
               θobs = Data[:, minimum(Array_iHeader): maximum(Array_iHeader)]

            # The depths were we have θ measurements
               Z = fill(0.0::Float64, Ndepth)

               i = 0::Int64
               for iHeader in Header
                  if occursin("Z=", iHeader)
                     i += 1
                     # Cleaning the header to get the integer
                     iHeader = replace(iHeader, "Z=" => "")
                     iHeader = replace(iHeader, "mm" => "")
                     iHeader = replace(iHeader, " " => "")
                     iHeader=  parse(Float64, iHeader)
                     Z[i] = iHeader
                  end # occursin("Z=", iHeader)
               end #  iHeader

         # REDUCING THE NUMBER OF SIMULATIONS SUCH THAT IT IS WITHIN THE SELECTED RANGE
      
         # # Compared to observed climate data
            Date_SimStart_θobs = DateTime(Year[1], Month[1], Day[1], Hour[1], Minute[1], Second[1])
            Date_End_θobs = DateTime(Year[end], Month[end], Day[end], Hour[end], Minute[end], Second[end])

            dateHypix.Date_SimStart = max(dateHypix.Date_SimStart, Date_SimStart_θobs)

         # SELECTING THE DATA WITHING FEASIBLE RANGE
            True = falses(Nit) # Initiating with false
            Date = fill(dateHypix.Date_SimStart::DateTime, Nit)
            iCount = 0 ::Int64
            for iT=1:Nit
               Date[iT] = DateTime(Year[iT], Month[iT], Day[iT], Hour[iT], Minute[iT], Second[iT])

               if (dateHypix.Date_SimStart ≤ Date[iT] ≤ dateHypix.Date_End)
                  iCount += 1
                  True[iT] = true
               end  # if
            end # for iT=1:Nit

            # New reduced number of simulations selected with dates
               Date = Date[True[1:Nit]]
               θobs = θobs[True[1:Nit], 1:Ndepth]

               Nit = count(True[1:Nit]) # New number of data

               if Nit == 0
                  error("\n ***** HyPix θDATA Date has no data *** \n")
               end

            # Updating
               dateHypix.Date_SimStart = Date[1]

         # REDUCING THE AMOUNT OF DATA TO HOURLY TODO

         # ∑T_Reduced = collect(range(Date[1], step=Hour[1], stop=Date[end])) 
            # ΔTimeStep = param.hyPix.ΔT_Output
            # if optionHypix.θobs_Reduced && ΔTimeStep < 86400
            # 	True = falses(Nit)
            # 	iCount = 0 
            # 	for iT=1:Nit
            # 		if hour(Date[iT]) == 0 && minute(Date[iT]) == 0
            # 			True[iT] = true
            # 			iCount += 1
            # 		end # if
            # 	end # for
            
            # 	# New reduced number of simulations selected with dates
            # 	Date = Date[True[1:Nit]]
            # 	θobs = θobs[True[1:Nit],1:Ndepth]
            # 	Nit = count(True[1:Nit]) # New reduced amount of data
            # end # θobs_Reduced

         # This will be computed at PrioProcess
            ∑T        = fill(0.0::Float64, Nit)
            ithetaObs = fill(0::Int64, Ndepth)

         # SAVING SPACE 
            # Data = nothing
            # True = nothing

         # STRUCTURE
   return dateHypix, θOBSERVATION(Date, Z, ithetaObs, Nit, Ndepth, θobs, ∑T) 
   end  # function: TIME_SERIES
   # ------------------------------------------------------------------
   
end  # module: readHypix
# ............................................................