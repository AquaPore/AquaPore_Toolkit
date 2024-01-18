# =============================================================
#		module: including
# =============================================================
# module including

using Suppressor

   @suppress begin
      include("Option.jl")
      include("Path.jl")
      include("Tool.jl")
      include("Cst.jl")
      include("Param.jl")

      include("Hydro/ΨminΨmax.jl")
      include("Hydro/HydroStruct.jl")
      include("HyPix/HorizonLayer.jl")
      include("Hydro/HydroRelation.jl")
      include("Hydro/Wrc.jl")
      include("Hydro/Kunsat.jl")
      include("Stats.jl")
      
      include("Psd/Psd_θr.jl")
      include("Psd/Psd_START.jl")

      include("Optim/Optimize.jl")
      include("Table.jl")
      include("Ksmodel/Struct_Ksmodel.jl")
      include("KsModel/θψ_2_KsModel.jl")
      include("Reading.jl")
      include("Distribution.jl")

      include("Sorptivity/Sorptivity.jl")
      include("Infiltration/BestFunc.jl")            
      include("Infiltration/Infiltration_START.jl")
      
      include("Plot.jl")

      include("Ksmodel/Start_KsModel.jl")

      include("Checking.jl")
      include("RockFragment/RockFragment.jl")
      include("HyPix/Veg/VegStruct.jl")
      include("HyPix/Discretisation.jl")

      include("HydroLab/OfHydrolab.jl")
      include("HydroLab/OptIndivSoil.jl")
      include("HydroLab/OptAllSoil.jl")
      include("HydroLab/HydrolabOpt.jl")

      include("NoCore/Smap/ReadSmap.jl")
      include("NoCore/Smap/TableSmap.jl")
      include("NoCore/Smap/Smap2Hypix.jl")
      include("NoCore/Smap/PlotSmap.jl")
      

      # include("HyPix/HypixStart.jl")
      # include("Temporary/Ks_Smap.jl")
      # include("NoCore/NSDR/ReadNsdr.jl")
   end # Suppressor

# end