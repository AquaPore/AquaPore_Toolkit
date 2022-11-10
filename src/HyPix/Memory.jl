# =============================================================
#		module: memory

# =============================================================
module memory
   export MEMORY_SCENARIO, MEMORY
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : MEMORY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function MEMORY(clim, N_∑T_Climate::Int64, Nz::Int64, obsθ, optionHypix, paramHypix)

         N_Memory = ceil(Int, N_∑T_Climate / paramHypix.ΔT_Min) + 10
         
         ΔEvaporation = fill(0.0::Float64, N_Memory)
         Hpond        = fill(0.0::Float64, N_Memory)
         ΔPet         = fill(0.0::Float64, N_Memory)
         ΔPrThroughfall          = fill(0.0::Float64, N_Memory)
         ΔRunoff      = fill(0.0::Float64, N_Memory)
         ΔT           = fill(0.0::Float64, N_Memory)
         ∑Pet         = fill(0.0::Float64, N_Memory)
         ∑PrThroughfall          = fill(0.0::Float64, N_Memory)
         ∑T           = fill(0.0::Float64, N_Memory)

         ΔSink = fill(0.0::Float64, N_Memory, Nz)
         Ψ     = Array{Float64}(undef, N_Memory, Nz)
         θ     = Array{Float64}(undef, N_Memory, Nz)
         Q     = Array{Float64}(undef, N_Memory, Nz+1)
         
         Pkₐᵥₑᵣ   = fill(1.0::Float64, Nz)
         Residual = fill(0.0::Float64,  Nz)
         ΔLnΨmax  = fill(0.0::Float64,  Nz)
         Ψ_Min    = fill(0.0::Float64,  Nz)
         Ψbest    = fill(0.0::Float64,  Nz)
         ∂K∂Ψ     = fill(0.0::Float64, Nz)
         ∂R∂Ψ     = fill(0.0::Float64, Nz)
         ∂R∂Ψ△    = fill(0.0::Float64, Nz)
         ∂R∂Ψ▽    = fill(0.0::Float64, Nz)

         K_Aver_Vect  = fill(1.0::Float64, Nz+1)
         K_Aver₀_Vect = fill(1.0::Float64, Nz+1)

         # Laiᵀ= fill(0.0::Float64, clim.N_Climate)
         CropCoeficientᵀ = fill(0.0::Float64, clim.N_Climate)

         if optionHypix.θobs
            θSim = fill(0.0::Float64, obsθ.Nit, Nz)
         else
            θSim = fill(0.0::Float64, 1, Nz)
         end
            
      return ∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑PrThroughfall, ∑T, CropCoeficientᵀ, Hpond, K_Aver_Vect, K_Aver₀_Vect, Pkₐᵥₑᵣ, Q, Residual, ΔEvaporation, ΔLnΨmax, ΔPet, ΔPrThroughfall, ΔRunoff, ΔSink, ΔT, θ, θSim, Ψ, Ψ_Min, Ψbest
      end  # function: MEMORY
   #.................................................................


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : MEMORY_SCENARIO
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function MEMORY_SCENARIO(N_Scenario, optionHypix, paramHypix)
         if optionHypix.opt.Optimisation
            # Nit_Reduced = (paramHypix.opt.iOptMultiStep_End - paramHypix.opt.iOptMultiStep_Start + 1) * N_Scenario

            Nit_Reduced = (paramHypix.opt.iOptMultiStep_End - paramHypix.opt.iOptMultiStep_Start + 1)
         else
            Nit_Reduced = N_Scenario
         end
        
        println("Nit_Reduced = $Nit_Reduced")
        println("N_Scenario = $N_Scenario")

         CccBest                    = fill(0.0::Float64, Nit_Reduced)
         Efficiency                 = fill(0.0::Float64, Nit_Reduced)
         Global_WaterBalance        = fill(0.0::Float64, Nit_Reduced)
         Global_WaterBalance_NormPr = fill(0.0::Float64, Nit_Reduced)
         NseBest                    = fill(0.0::Float64, Nit_Reduced)
         WilmotBest                 = fill(0.0::Float64, Nit_Reduced)
         WofBest                    = fill(0.0::Float64, Nit_Reduced)
         iNonConverge_iOpt          = fill(0.0::Float64, Nit_Reduced)
         ΔRunTimeHypix              = fill(0.0::Float64, Nit_Reduced)
         ΔT_Average                 = fill(0.0::Float64, Nit_Reduced)
         θroot_Mean                 = fill(0.0::Float64, Nit_Reduced)
         ∑Pet_Net                   = fill(0.0::Float64, Nit_Reduced)
         ∑Pr_Clim                   = fill(0.0::Float64, Nit_Reduced)
         ∑Q_Z                       = fill(0.0::Float64, Nit_Reduced)
         ∑ΔQ_Bot                    = fill(0.0::Float64, Nit_Reduced)
         ∑∑ΔSink                    = fill(0.0::Float64, Nit_Reduced)
         
      return ∑∑ΔSink, ∑Pet_Net, ∑Pr_Clim, ∑Q_Z, ∑ΔQ_Bot, CccBest, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iOpt, NseBest, WilmotBest, WofBest, ΔRunTimeHypix, ΔT_Average, θroot_Mean
      end  # function: MEMORY_STEOPT
   # ------------------------------------------------------------------

end  # module: memory 

# ............................................................