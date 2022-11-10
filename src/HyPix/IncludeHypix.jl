using Suppressor

Path_SoilWater = dirname(dirname(@__DIR__)) * "/src/" # moving down the path twice

# @suppress begin	
   include("../Cst.jl")
   include("../Hydro/Wrc.jl")
   include("../Hydro/Kunsat.jl")
   include("../Tool.jl")
   include("../Hydro/HydroStruct.jl")
   include("../Hydro/HydroRelation.jl")
   include("../Stats.jl")
   include("HorizonLayer.jl")

   include("ReadHypix/ReadLinkingFile.jl")
   include("ReadHypix/OptionsHypix.jl")
   include("ReadHypix/ParamsHypix.jl")
   include("ReadHypix/PathsHypix.jl")
   include("DatesHypix.jl")
   include("Opt/ΘObs.jl")
   include("Climate.jl")
   include("Discretisation.jl")
   include("Memory.jl")
   include("Veg/VegStruct.jl")
   include("HypixStats/θaver.jl")
   include("OutputHypix/TableHypix.jl")
   include("HydroSmooth.jl")

   include("Interpolate.jl")
   include("ReadHypix/ReadHypix.jl")
   include("WaterBalance.jl")
   include("Veg/Interception.jl")
   include("Veg/RootWaterUptake.jl")
   include("../Hydro/ΨminΨmax.jl")
   include("TimeStep.jl")
   include("Evaporation.jl")
   include("Opt/OfHypix.jl")
   include("Veg/Pet.jl")
   include("../Sorptivity/Sorptivity.jl")
   include("Ponding.jl")
   include("Flux.jl")
   include("Residual.jl")
   include("Richard.jl")
   include("HypixModel.jl")
   include("Opt/HypixOpt.jl")
   include("HypixStats/ΔΔtchange.jl")

   include("OutputHypix/PlotHypix.jl")
   # include("OutputHypix/PlotOther.jl")

# end # @suppress begin	
