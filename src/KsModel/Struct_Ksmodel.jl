# =============================================================
#		module: ksModel
# =============================================================
module ksModel

   Base.@kwdef mutable struct KSMODELτ
      τ₁ₐ             :: Vector{Float64}
      τ₁ₐ_Max         :: Vector{Float64}
      τ₁ₐ_Min         :: Vector{Float64}
         
      τ₁ₐMac          :: Vector{Float64}
      τ₁ₐMac_Max      :: Vector{Float64}
      τ₁ₐMac_Min      :: Vector{Float64}
      
      τ₂ₐ             :: Vector{Float64}
      τ₂ₐ_Max         :: Vector{Float64}
      τ₂ₐ_Min         :: Vector{Float64}
      
      τ₂ₐMac     :: Vector{Float64}
      τ₂ₐMac_Max :: Vector{Float64}
      τ₂ₐMac_Min :: Vector{Float64}
      
      τ₃ₐ             :: Vector{Float64}
      τ₃ₐ_Min         :: Vector{Float64}
      τ₃ₐ_Max         :: Vector{Float64}
   
      τ₃ₐMac          :: Vector{Float64}
      τ₃ₐMac_Max      :: Vector{Float64}
      τ₃ₐMac_Min      :: Vector{Float64}
         
      τclay₀          :: Vector{Float64}
      τclay₀_Max      :: Vector{Float64}
      τclay₀_Min      :: Vector{Float64}
      
      τclayₘₐₓ        :: Vector{Float64}
      τclayₘₐₓ_Max    :: Vector{Float64}
      τclayₘₐₓ_Min    :: Vector{Float64}

      τclayΔθsr       :: Vector{Float64}
      τclayΔθsr_Max   :: Vector{Float64}
      τclayΔθsr_Min   :: Vector{Float64}
         
      Nse_τ           :: Vector{Float64}
      Rmse_τ          :: Vector{Float64}
      Wilmot_τ        :: Vector{Float64}
      Ccc_τ           :: Vector{Float64}
   end # mutable struct KSMODEL


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : STRUCT_KSMODEL
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function STRUCT_KSMODEL(; N_KsClass=N_KsClass::Int64)

         τ₁ₐ           = fill(0.0::Float64,  N_KsClass)
         τ₁ₐ_Max       = fill(0.0::Float64,  N_KsClass)
         τ₁ₐ_Min       = fill(0.0::Float64,  N_KsClass)

         τ₁ₐMac        = fill(0.0::Float64,  N_KsClass)
         τ₁ₐMac_Max    = fill(0.0::Float64,  N_KsClass)
         τ₁ₐMac_Min    = fill(0.0::Float64,  N_KsClass)

         τ₂ₐ           = fill(0.0::Float64,   N_KsClass)
         τ₂ₐ_Max       = fill(0.0::Float64,   N_KsClass)
         τ₂ₐ_Min       = fill(0.0::Float64,   N_KsClass)

         τ₂ₐMac        = fill(0.0::Float64,   N_KsClass)
         τ₂ₐMac_Max    = fill(0.0::Float64,  N_KsClass)
         τ₂ₐMac_Min    = fill(0.0::Float64,  N_KsClass)
         
         τ₃ₐ           = fill(0.0::Float64,   N_KsClass)
         τ₃ₐ_Max       = fill(0.0::Float64,  N_KsClass)
         τ₃ₐ_Min       = fill(0.0::Float64,  N_KsClass)
         
         τ₃ₐMac        = fill(0.0::Float64,  N_KsClass)
         τ₃ₐMac_Max    = fill(0.0::Float64,  N_KsClass)
         τ₃ₐMac_Min    = fill(0.0::Float64,  N_KsClass)

         τclay₀        = fill(0.0::Float64,  N_KsClass)
         τclay₀_Max    = fill(0.0::Float64,  N_KsClass)
         τclay₀_Min    = fill(0.0::Float64,  N_KsClass)

         τclayₘₐₓ      = fill(0.0::Float64,   N_KsClass)
         τclayₘₐₓ_Max  = fill(0.0::Float64,  N_KsClass)
         τclayₘₐₓ_Min  = fill(0.0::Float64,  N_KsClass)

         τclayΔθsr     = fill(0.0::Float64,   N_KsClass)
         τclayΔθsr_Max = fill(0.0::Float64,   N_KsClass)
         τclayΔθsr_Min = fill(0.0::Float64,   N_KsClass)

         Nse_τ         = fill(0.0::Float64,   N_KsClass)
         Rmse_τ        = fill(0.0::Float64,   N_KsClass)
         Wilmot_τ      = fill(0.0::Float64,   N_KsClass)
         Ccc_τ         = fill(0.0::Float64,   N_KsClass)

      return ksmodelτ = KSMODELτ(
         τ₁ₐ,
         τ₁ₐ_Max,
         τ₁ₐ_Min,
         τ₁ₐMac, 
         τ₁ₐMac_Max,
         τ₁ₐMac_Min,
         τ₂ₐ,
         τ₂ₐ_Max,
         τ₂ₐ_Min,
         τ₂ₐMac,
         τ₂ₐMac_Max,
         τ₂ₐMac_Min,
         τ₃ₐ,
         τ₃ₐ_Min,
         τ₃ₐ_Max,
         τ₃ₐMac,
         τ₃ₐMac_Max,
         τ₃ₐMac_Min,
         τclay₀,
         τclay₀_Max,
         τclay₀_Min,
         τclayₘₐₓ,
         τclayₘₐₓ_Max,
         τclayₘₐₓ_Min,
         τclayΔθsr,
         τclayΔθsr_Max,
         τclayΔθsr_Min,
         Nse_τ,
         Rmse_τ,
         Wilmot_τ,
         Ccc_τ,
      )

      end  # function: STRUCT_KSMODEL

   end  # module: ksModel
# ............................................................
