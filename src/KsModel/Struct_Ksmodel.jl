# =============================================================
#		module: ksModel
# =============================================================
module ksModel

   Base.@kwdef mutable struct KSMODELτ
      τ₁ₐ::Vector{Float64}
      τ₁ᵦ::Vector{Float64}
      τ₂ₐ::Vector{Float64}
      τ₂ᵦ::Vector{Float64}
      τ₃ₐ::Vector{Float64}
      τ₃ᵦ::Vector{Float64}
      τ₁ₐMac::Vector{Float64}
      τ₁ᵦMac::Vector{Float64}
      τ₂ₐMac::Vector{Float64}
      τ₂ᵦMac::Vector{Float64}
      τ₃ₐMac::Vector{Float64}
      τ₃ᵦMac::Vector{Float64}

      τ₁ₐ_Min::Vector{Float64}
      τ₁ᵦ_Min::Vector{Float64}
      τ₂ₐ_Min::Vector{Float64}
      τ₂ᵦ_Min::Vector{Float64}
      τ₃ₐ_Min::Vector{Float64}
      τ₃ᵦ_Min::Vector{Float64}
      τ₁ₐMac_Min::Vector{Float64}
      τ₁ᵦMac_Min::Vector{Float64}
      τ₂ₐMac_Min::Vector{Float64}
      τ₂ᵦMac_Min::Vector{Float64}
      τ₃ₐMac_Min::Vector{Float64}
      τ₃ᵦMac_Min::Vector{Float64}

      τ₁ₐ_Max::Vector{Float64}
      τ₁ᵦ_Max::Vector{Float64}
      τ₂ₐ_Max::Vector{Float64}
      τ₂ᵦ_Max::Vector{Float64}
      τ₃ₐ_Max::Vector{Float64}
      τ₃ᵦ_Max::Vector{Float64}
      τ₁ₐMac_Max::Vector{Float64}
      τ₁ᵦMac_Max::Vector{Float64}
      τ₂ₐMac_Max::Vector{Float64}
      τ₂ᵦMac_Max::Vector{Float64}
      τ₃ₐMac_Max::Vector{Float64}
      τ₃ᵦMac_Max::Vector{Float64}

      Nse_τ::Vector{Float64}
      Rmse_τ::Vector{Float64}
      Wilmot_τ::Vector{Float64}
      Ccc_τ::Vector{Float64}
   end # mutable struct KSMODEL


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : STRUCT_KSMODEL
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function STRUCT_KSMODEL(; Nτ_Layer=Nτ_Layer::Int64)

         τ₁ₐ        = fill(0.0::Float64, Nτ_Layer)
         τ₁ᵦ        = fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ        = fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ        = fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ        = fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ        = fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac     = fill(0.0::Float64, Nτ_Layer)
         τ₁ᵦMac     = fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac     = fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac     = fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac     = fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac     = fill(0.0::Float64, Nτ_Layer)

         τ₁ₐ_Min    = fill(0.0::Float64, Nτ_Layer)
         τ₁ᵦ_Min    = fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_Min    = fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_Min    = fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_Min    = fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_Min    = fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_Min = fill(0.0::Float64, Nτ_Layer)
         τ₁ᵦMac_Min = fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_Min = fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_Min = fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_Min = fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_Min = fill(0.0::Float64, Nτ_Layer)

         τ₁ₐ_Max    = fill(0.0::Float64, Nτ_Layer)
         τ₁ᵦ_Max    = fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_Max    = fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_Max    = fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_Max    = fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_Max    = fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_Max = fill(0.0::Float64, Nτ_Layer)
         τ₁ᵦMac_Max = fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_Max = fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_Max = fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_Max = fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_Max = fill(0.0::Float64, Nτ_Layer)

         Nse_τ      = fill(0.0::Float64, Nτ_Layer)
         Rmse_τ     = fill(0.0::Float64, Nτ_Layer)
         Wilmot_τ   = fill(0.0::Float64, Nτ_Layer)
         Ccc_τ      = fill(0.0::Float64, Nτ_Layer)

      return ksmodelτ = KSMODELτ(τ₁ₐ,τ₁ᵦ,τ₂ₐ,τ₂ᵦ,τ₃ₐ,τ₃ᵦ,τ₁ₐMac,τ₁ᵦMac,τ₂ₐMac,τ₂ᵦMac,τ₃ₐMac,τ₃ᵦMac,τ₁ₐ_Min,τ₁ᵦ_Min,τ₂ₐ_Min,τ₂ᵦ_Min,τ₃ₐ_Min,τ₃ᵦ_Min,τ₁ₐMac_Min,τ₁ᵦMac_Min,τ₂ₐMac_Min,τ₂ᵦMac_Min,τ₃ₐMac_Min,τ₃ᵦMac_Min,τ₁ₐ_Max,τ₁ᵦ_Max,τ₂ₐ_Max,τ₂ᵦ_Max,τ₃ₐ_Max,τ₃ᵦ_Max,τ₁ₐMac_Max,τ₁ᵦMac_Max,τ₂ₐMac_Max,τ₂ᵦMac_Max,τ₃ₐMac_Max,τ₃ᵦMac_Max, Nse_τ,Rmse_τ,Wilmot_τ,Ccc_τ)
      end  # function: STRUCT_KSMODEL

   end  # module: ksModel
# ............................................................
