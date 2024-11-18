# =============================================================
#		module: ksModel
# =============================================================
module ksModel

   Base.@kwdef mutable struct KSMODELτ

      τ₁ₐ             :: Vector{Float64}
      τ₁ₐ_Max         :: Vector{Float64}
      τ₁ₐ_Min         :: Vector{Float64}

      τ₁ₐ_1           :: Vector{Float64}
      τ₁ₐ_1_Max       :: Vector{Float64}
      τ₁ₐ_1_Min       :: Vector{Float64}
      
      τ₁ₐ_2           :: Vector{Float64}
      τ₁ₐ_2_Max       :: Vector{Float64}
      τ₁ₐ_2_Min       :: Vector{Float64}

      τ₁ₐ_3           :: Vector{Float64}
      τ₁ₐ_3_Max       :: Vector{Float64}
      τ₁ₐ_3_Min       :: Vector{Float64}
      
      τ₁ₐ_4           :: Vector{Float64}
      τ₁ₐ_4_Max       :: Vector{Float64}
      τ₁ₐ_4_Min       :: Vector{Float64}
      
      τ₁ₐ_5           :: Vector{Float64}
      τ₁ₐ_5_Max       :: Vector{Float64}
      τ₁ₐ_5_Min       :: Vector{Float64}
            
      τ₁ₐMac          :: Vector{Float64}
      τ₁ₐMac_Max      :: Vector{Float64}
      τ₁ₐMac_Min      :: Vector{Float64}
      
      τ₁ₐMac_1        :: Vector{Float64}
      τ₁ₐMac_1_Max    :: Vector{Float64}
      τ₁ₐMac_1_Min    :: Vector{Float64}
      
      τ₁ₐMac_2        :: Vector{Float64}
      τ₁ₐMac_2_Max    :: Vector{Float64}
      τ₁ₐMac_2_Min    :: Vector{Float64}
      
      τ₁ₐMac_3        :: Vector{Float64}
      τ₁ₐMac_3_Max    :: Vector{Float64}
      τ₁ₐMac_3_Min    :: Vector{Float64}
      
      τ₁ₐMac_4        :: Vector{Float64}
      τ₁ₐMac_4_Max    :: Vector{Float64}
      τ₁ₐMac_4_Min    :: Vector{Float64}
      
      τ₁ₐMac_5        :: Vector{Float64}
      τ₁ₐMac_5_Max    :: Vector{Float64}
      τ₁ₐMac_5_Min    :: Vector{Float64}
            
      τ₂ₐ             :: Vector{Float64}
      τ₂ₐ_Max         :: Vector{Float64}
      τ₂ₐ_Min         :: Vector{Float64}
      
      τ₂ₐ_1           :: Vector{Float64}
      τ₂ₐ_1_Max       :: Vector{Float64}
      τ₂ₐ_1_Min       :: Vector{Float64}
      
      τ₂ₐ_2           :: Vector{Float64}
      τ₂ₐ_2_Max       :: Vector{Float64}
      τ₂ₐ_2_Min       :: Vector{Float64}
      
      τ₂ₐ_3           :: Vector{Float64}
      τ₂ₐ_3_Max       :: Vector{Float64}
      τ₂ₐ_3_Min       :: Vector{Float64}
      
      τ₂ₐ_4           :: Vector{Float64}
      τ₂ₐ_4_Max       :: Vector{Float64}
      τ₂ₐ_4_Min       :: Vector{Float64}
      
      τ₂ₐ_5      :: Vector{Float64}
      τ₂ₐ_5_Max  :: Vector{Float64}
      τ₂ₐ_5_Min  :: Vector{Float64}
      
      τ₂ₐMac     :: Vector{Float64}
      τ₂ₐMac_Max :: Vector{Float64}
      τ₂ₐMac_MIn :: Vector{Float64}
      
      τ₂ₐMac_1        :: Vector{Float64}
      τ₂ₐMac_1_Max    :: Vector{Float64}
      τ₂ₐMac_1_Min    :: Vector{Float64}
      
      τ₂ₐMac_2        :: Vector{Float64}
      τ₂ₐMac_2_Max    :: Vector{Float64}
      τ₂ₐMac_2_Min    :: Vector{Float64}
      
      τ₂ₐMac_3        :: Vector{Float64}
      τ₂ₐMac_3_Max    :: Vector{Float64}
      τ₂ₐMac_3_Min    :: Vector{Float64}
      
      τ₂ₐMac_4        :: Vector{Float64}
      τ₂ₐMac_4_Max    :: Vector{Float64}
      τ₂ₐMac_4_Min    :: Vector{Float64}
      
      τ₂ₐMac_5        :: Vector{Float64}
      τ₂ₐMac_5_Max    :: Vector{Float64}
      τ₂ₐMac_5_Min    :: Vector{Float64}
      τ₂ₐMac_Min      :: Vector{Float64}
      
      τ₃ₐ             :: Vector{Float64}
      τ₃ₐ_Min         :: Vector{Float64}
      τ₃ₐ_Max         :: Vector{Float64}
      
      τ₃ₐ_1           :: Vector{Float64}
      τ₃ₐ_1_Max       :: Vector{Float64}
      τ₃ₐ_1_Min       :: Vector{Float64}
      
      τ₃ₐ_2           :: Vector{Float64}
      τ₃ₐ_2_Max       :: Vector{Float64}
      τ₃ₐ_2_Min       :: Vector{Float64}
      
      τ₃ₐ_3           :: Vector{Float64}
      τ₃ₐ_3_Max       :: Vector{Float64}
      τ₃ₐ_3_Min       :: Vector{Float64}
      
      τ₃ₐ_4           :: Vector{Float64}
      τ₃ₐ_4_Max       :: Vector{Float64}
      τ₃ₐ_4_Min       :: Vector{Float64}
      
      τ₃ₐ_5           :: Vector{Float64}
      τ₃ₐ_5_Max       :: Vector{Float64}
      τ₃ₐ_5_Min       :: Vector{Float64}
      
      τ₃ₐMac          :: Vector{Float64}
      τ₃ₐMac_Max      :: Vector{Float64}
      τ₃ₐMac_Min      :: Vector{Float64}
      
      τ₃ₐMac_1        :: Vector{Float64}
      τ₃ₐMac_1_Max    :: Vector{Float64}
      τ₃ₐMac_1_Min    :: Vector{Float64}
      
      τ₃ₐMac_2        :: Vector{Float64}
      τ₃ₐMac_2_Max    :: Vector{Float64}
      τ₃ₐMac_2_Min    :: Vector{Float64}
      
      τ₃ₐMac_3        :: Vector{Float64}
      τ₃ₐMac_3_Max    :: Vector{Float64}
      τ₃ₐMac_3_Min    :: Vector{Float64}
      
      τ₃ₐMac_4        :: Vector{Float64}
      τ₃ₐMac_4_Max    :: Vector{Float64}
      τ₃ₐMac_4_Min    :: Vector{Float64}
      
      τ₃ₐMac_5        :: Vector{Float64}
      τ₃ₐMac_5_Max    :: Vector{Float64}
      τ₃ₐMac_5_Min    :: Vector{Float64}
         
      τclay₀          :: Vector{Float64}
      τclay₀_Max      :: Vector{Float64}
      τclay₀_Min      :: Vector{Float64}

      τclay₀_1        :: Vector{Float64}
      τclay₀_1_Max    :: Vector{Float64}
      τclay₀_1_Min    :: Vector{Float64}
      
      τclay₀_2        :: Vector{Float64}
      τclay₀_2_Max    :: Vector{Float64}
      τclay₀_2_Min    :: Vector{Float64}
      
      τclay₀_3        :: Vector{Float64}
      τclay₀_3_Max    :: Vector{Float64}
      τclay₀_3_Min    :: Vector{Float64}
      
      τclay₀_4        :: Vector{Float64}
      τclay₀_4_Max    :: Vector{Float64}
      τclay₀_4_Min    :: Vector{Float64}
      
      τclay₀_5        :: Vector{Float64}
      τclay₀_5_Max    :: Vector{Float64}
      τclay₀_5_Min    :: Vector{Float64}
      
      τclayₘₐₓ        :: Vector{Float64}
      τclayₘₐₓ_Max    :: Vector{Float64}
      τclayₘₐₓ_Min    :: Vector{Float64}

      τclayₘₐₓ_1      :: Vector{Float64}
      τclayₘₐₓ_1_Max  :: Vector{Float64}
      τclayₘₐₓ_1_Min  :: Vector{Float64}
   
      τclayₘₐₓ_2      :: Vector{Float64}
      τclayₘₐₓ_2_Max  :: Vector{Float64}
      τclayₘₐₓ_2_Min  :: Vector{Float64}
   
      τclayₘₐₓ_3      :: Vector{Float64}
      τclayₘₐₓ_3_Max  :: Vector{Float64}
      τclayₘₐₓ_3_Min  :: Vector{Float64}
      
      τclayₘₐₓ_4      :: Vector{Float64}
      τclayₘₐₓ_4_Max  :: Vector{Float64}
      τclayₘₐₓ_4_Min  :: Vector{Float64}

      τclayₘₐₓ_5      :: Vector{Float64}
      τclayₘₐₓ_5_Max  :: Vector{Float64}
      τclayₘₐₓ_5_Min  :: Vector{Float64}

      τclayΔθsr       :: Vector{Float64}
      τclayΔθsr_Max   :: Vector{Float64}
      τclayΔθsr_Min   :: Vector{Float64}

      τclayΔθsr_1     :: Vector{Float64}
      τclayΔθsr_1_Max :: Vector{Float64}
      τclayΔθsr_1_Min :: Vector{Float64}

      τclayΔθsr_2     :: Vector{Float64}
      τclayΔθsr_2_Max :: Vector{Float64}
      τclayΔθsr_2_Min :: Vector{Float64}

      τclayΔθsr_3     :: Vector{Float64}
      τclayΔθsr_3_Max :: Vector{Float64}
      τclayΔθsr_3_Min :: Vector{Float64}

      τclayΔθsr_4     :: Vector{Float64}
      τclayΔθsr_4_Max :: Vector{Float64}
      τclayΔθsr_4_Min :: Vector{Float64}

      τclayΔθsr_5     :: Vector{Float64}
      τclayΔθsr_5_Max :: Vector{Float64}
      τclayΔθsr_5_Min :: Vector{Float64}
         
      Nse_τ           :: Vector{Float64}
      Rmse_τ          :: Vector{Float64}
      Wilmot_τ        :: Vector{Float64}
      Ccc_τ           :: Vector{Float64}
   end # mutable struct KSMODEL


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : STRUCT_KSMODEL
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function STRUCT_KSMODEL(; Nτ_Layer=Nτ_Layer::Int64)

         τ₁ₐ                = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_1              = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_1_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_1_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_2              = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_2_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_2_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_3              = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_3_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_3_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_4              = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_4_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_4_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_5              = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_5_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_5_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_Max            = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐ_Min            = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac             = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_1           = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_1_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_1_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_2           = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_2_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_2_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_3           = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_3_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_3_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_4           = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_4_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_4_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_5           = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_5_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_5_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_Max         = fill(0.0::Float64,  Nτ_Layer)
         τ₁ₐMac_Min         = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ                = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_1              = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_1_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_1_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_2              = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_2_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_2_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_3              = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_3_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_3_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_4              = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_4_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_4_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_5              = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_5_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_5_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_Max            = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐ_Min            = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac             = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_1           = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_1_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_1_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_2           = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_2_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_2_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_3           = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_3_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_3_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_4           = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_4_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_4_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_5           = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_5_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_5_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_Max         = fill(0.0::Float64,  Nτ_Layer)
         τ₂ₐMac_Min         = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ                = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_1              = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_1_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_1_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_2              = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_2_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_2_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_3              = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_3_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_3_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_4              = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_4_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_4_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_5              = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_5_Max          = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_5_Min          = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_Max            = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐ_Min            = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac             = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_1           = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_1_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_1_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_2           = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_2_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_2_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_3           = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_3_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_3_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_4           = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_4_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_4_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_5           = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_5_Max       = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_5_Min       = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_Max         = fill(0.0::Float64,  Nτ_Layer)
         τ₃ₐMac_Min         = fill(0.0::Float64,  Nτ_Layer)
         τclay₀             = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_1           = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_1_Max       = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_1_Min       = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_2           = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_2_Max       = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_2_Min       = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_3           = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_3_Max       = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_3_Min       = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_4           = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_4_Max       = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_4_Min       = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_5           = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_5_Max       = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_5_Min       = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_Max         = fill(0.0::Float64,  Nτ_Layer)
         τclay₀_Min         = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ           = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_1         = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_1_Max     = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_1_Min     = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_2         = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_2_Max     = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_2_Min     = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_3         = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_3_Max     = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_3_Min     = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_4         = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_4_Max     = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_4_Min     = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_5         = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_5_Max     = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_5_Min     = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_Max       = fill(0.0::Float64,  Nτ_Layer)
         τclayₘₐₓ_Min       = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr          = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_1        = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_1_Max    = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_1_Min    = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_2        = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_2_Max    = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_2_Min    = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_3        = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_3_Max    = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_3_Min    = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_4        = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_4_Max    = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_4_Min    = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_5        = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_5_Max    = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_5_Min    = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_Max      = fill(0.0::Float64,  Nτ_Layer)
         τclayΔθsr_Min      = fill(0.0::Float64,  Nτ_Layer)

         Nse_τ      = fill(0.0::Float64,  Nτ_Layer)
         Rmse_τ     = fill(0.0::Float64,  Nτ_Layer)
         Wilmot_τ   = fill(0.0::Float64,  Nτ_Layer)
         Ccc_τ      = fill(0.0::Float64,  Nτ_Layer)

      return ksmodelτ = KSMODELτ(τ₁ₐ,
         τ₁ₐ_Max,
         τ₁ₐ_Min,
         τ₁ₐ_1, 
         τ₁ₐ_1_Max,
         τ₁ₐ_1_Min,
         τ₁ₐ_2, 
         τ₁ₐ_2_Max,
         τ₁ₐ_2_Min,
         τ₁ₐ_3, 
         τ₁ₐ_3_Max,
         τ₁ₐ_3_Min,
         τ₁ₐ_4, 
         τ₁ₐ_4_Max,
         τ₁ₐ_4_Min,
         τ₁ₐ_5, 
         τ₁ₐ_5_Max,
         τ₁ₐ_5_Min,
         τ₁ₐMac, 
         τ₁ₐMac_Max,
         τ₁ₐMac_Min,
         τ₁ₐMac_1,
         τ₁ₐMac_1_Max,
         τ₁ₐMac_1_Min,
         τ₁ₐMac_2,
         τ₁ₐMac_2_Max,
         τ₁ₐMac_2_Min,
         τ₁ₐMac_3,
         τ₁ₐMac_3_Max,
         τ₁ₐMac_3_Min,
         τ₁ₐMac_4,
         τ₁ₐMac_4_Max,
         τ₁ₐMac_4_Min,
         τ₁ₐMac_5,
         τ₁ₐMac_5_Max,
         τ₁ₐMac_5_Min,
         τ₂ₐ,
         τ₂ₐ_Max,
         τ₂ₐ_Min,
         τ₂ₐ_1,
         τ₂ₐ_1_Max,
         τ₂ₐ_1_Min,
         τ₂ₐ_2,
         τ₂ₐ_2_Max,
         τ₂ₐ_2_Min,
         τ₂ₐ_3,
         τ₂ₐ_3_Max,
         τ₂ₐ_3_Min,
         τ₂ₐ_4,
         τ₂ₐ_4_Max,
         τ₂ₐ_4_Min,
         τ₂ₐ_5,
         τ₂ₐ_5_Max,
         τ₂ₐ_5_Min,
         τ₂ₐMac,
         τ₂ₐMac_Max,
         τ₂ₐMac_Min,
         τ₂ₐMac_1,
         τ₂ₐMac_1_Max,
         τ₂ₐMac_1_Min,
         τ₂ₐMac_2,
         τ₂ₐMac_2_Max,
         τ₂ₐMac_2_Min,
         τ₂ₐMac_3,
         τ₂ₐMac_3_Max,
         τ₂ₐMac_3_Min,
         τ₂ₐMac_4,
         τ₂ₐMac_4_Max,
         τ₂ₐMac_4_Min,
         τ₂ₐMac_5,
         τ₂ₐMac_5_Max,
         τ₂ₐMac_5_Min,
         τ₂ₐMac_Min,
         τ₃ₐ,
         τ₃ₐ_Min,
         τ₃ₐ_Max,
         τ₃ₐ_1,
         τ₃ₐ_1_Max,
         τ₃ₐ_1_Min,
         τ₃ₐ_2,
         τ₃ₐ_2_Max,
         τ₃ₐ_2_Min,
         τ₃ₐ_3,
         τ₃ₐ_3_Max,
         τ₃ₐ_3_Min,
         τ₃ₐ_4,
         τ₃ₐ_4_Max,
         τ₃ₐ_4_Min,
         τ₃ₐ_5,
         τ₃ₐ_5_Max,
         τ₃ₐ_5_Min,
         τ₃ₐMac,
         τ₃ₐMac_Max,
         τ₃ₐMac_Min,
         τ₃ₐMac_1,
         τ₃ₐMac_1_Max,
         τ₃ₐMac_1_Min,
         τ₃ₐMac_2,
         τ₃ₐMac_2_Max,
         τ₃ₐMac_2_Min,
         τ₃ₐMac_3,
         τ₃ₐMac_3_Max,
         τ₃ₐMac_3_Min,
         τ₃ₐMac_4,
         τ₃ₐMac_4_Max,
         τ₃ₐMac_4_Min,
         τ₃ₐMac_5,
         τ₃ₐMac_5_Max,
         τ₃ₐMac_5_Min,
         τclay₀,
         τclay₀_Max,
         τclay₀_Min,
         τclay₀_1,
         τclay₀_1_Max,
         τclay₀_1_Min,
         τclay₀_2,
         τclay₀_2_Max,
         τclay₀_2_Min,
         τclay₀_3,
         τclay₀_3_Max,
         τclay₀_3_Min,
         τclay₀_4,
         τclay₀_4_Max,
         τclay₀_4_Min,
         τclay₀_5,
         τclay₀_5_Max,
         τclay₀_5_Min,
         τclayₘₐₓ,
         τclayₘₐₓ_Max,
         τclayₘₐₓ_Min,
         τclayₘₐₓ_1,
         τclayₘₐₓ_1_Max,
         τclayₘₐₓ_1_Min,
         τclayₘₐₓ_2,
         τclayₘₐₓ_2_Max,
         τclayₘₐₓ_2_Min,
         τclayₘₐₓ_3,
         τclayₘₐₓ_3_Max,
         τclayₘₐₓ_3_Min,
         τclayₘₐₓ_4,
         τclayₘₐₓ_4_Max,
         τclayₘₐₓ_4_Min,
         τclayₘₐₓ_5,
         τclayₘₐₓ_5_Max,
         τclayₘₐₓ_5_Min,
         τclayΔθsr,
         τclayΔθsr_Max,
         τclayΔθsr_Min,
         τclayΔθsr_1,
         τclayΔθsr_1_Max,
         τclayΔθsr_1_Min,
         τclayΔθsr_2,
         τclayΔθsr_2_Max,
         τclayΔθsr_2_Min,
         τclayΔθsr_3,
         τclayΔθsr_3_Max,
         τclayΔθsr_3_Min,
         τclayΔθsr_4,
         τclayΔθsr_4_Max,
         τclayΔθsr_4_Min,
         τclayΔθsr_5,
         τclayΔθsr_5_Max,
         τclayΔθsr_5_Min,
         Nse_τ,
         Rmse_τ,
         Wilmot_τ,
         Ccc_τ,
      )

      end  # function: STRUCT_KSMODEL

   end  # module: ksModel
# ............................................................
