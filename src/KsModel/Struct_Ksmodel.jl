# =============================================================
#		module: ksModel
# =============================================================
module ksModel

   Base.@kwdef mutable struct KSMODELτ
      τ₁ₐ::Vector{Float64}
      τclay₀::Vector{Float64}
      τ₂ₐ::Vector{Float64}
      τclayₘₐₓ::Vector{Float64}
      τ₃ₐ::Vector{Float64}
      τclayΔθsr::Vector{Float64}
      τ₁ₐMac::Vector{Float64}
      τclay₀Mac::Vector{Float64}
      τ₂ₐMac::Vector{Float64}
      τclayₘₐₓMac::Vector{Float64}
      τ₃ₐMac::Vector{Float64}
      τclayΔθsrMac::Vector{Float64}

      τ₁ₐ_1::Vector{Float64}
      τclay₀_1::Vector{Float64}
      τ₂ₐ_1::Vector{Float64}
      τclayₘₐₓ_1::Vector{Float64}
      τ₃ₐ_1::Vector{Float64}
      τclayΔθsr_1::Vector{Float64}
      τ₁ₐMac_1::Vector{Float64}
      τclay₀Mac_1::Vector{Float64}
      τ₂ₐMac_1::Vector{Float64}
      τclayₘₐₓMac_1::Vector{Float64}
      τ₃ₐMac_1::Vector{Float64}
      τclayΔθsrMac_1::Vector{Float64}

      τ₁ₐ_2::Vector{Float64}
      τclay₀_2::Vector{Float64}
      τ₂ₐ_2::Vector{Float64}
      τclayₘₐₓ_2::Vector{Float64}
      τ₃ₐ_2::Vector{Float64}
      τclayΔθsr_2::Vector{Float64}
      τ₁ₐMac_2::Vector{Float64}
      τclay₀Mac_2::Vector{Float64}
      τ₂ₐMac_2::Vector{Float64}
      τclayₘₐₓMac_2::Vector{Float64}
      τ₃ₐMac_2::Vector{Float64}
      τclayΔθsrMac_2::Vector{Float64}

      τ₁ₐ_3::Vector{Float64}
      τclay₀_3::Vector{Float64}
      τ₂ₐ_3::Vector{Float64}
      τclayₘₐₓ_3::Vector{Float64}
      τ₃ₐ_3::Vector{Float64}
      τclayΔθsr_3::Vector{Float64}
      τ₁ₐMac_3::Vector{Float64}
      τclay₀Mac_3::Vector{Float64}
      τ₂ₐMac_3::Vector{Float64}
      τclayₘₐₓMac_3::Vector{Float64}
      τ₃ₐMac_3::Vector{Float64}
      τclayΔθsrMac_3::Vector{Float64}

      τ₁ₐ_4::Vector{Float64}
      τclay₀_4::Vector{Float64}
      τ₂ₐ_4::Vector{Float64}
      τclayₘₐₓ_4::Vector{Float64}
      τ₃ₐ_4::Vector{Float64}
      τclayΔθsr_4::Vector{Float64}
      τ₁ₐMac_4::Vector{Float64}
      τclay₀Mac_4::Vector{Float64}
      τ₂ₐMac_4::Vector{Float64}
      τclayₘₐₓMac_4::Vector{Float64}
      τ₃ₐMac_4::Vector{Float64}
      τclayΔθsrMac_4::Vector{Float64}

      τ₁ₐ_5::Vector{Float64}
      τclay₀_5::Vector{Float64}
      τ₂ₐ_5::Vector{Float64}
      τclayₘₐₓ_5::Vector{Float64}
      τ₃ₐ_5::Vector{Float64}
      τclayΔθsr_5::Vector{Float64}
      τ₁ₐMac_5::Vector{Float64}
      τclay₀Mac_5::Vector{Float64}
      τ₂ₐMac_5::Vector{Float64}
      τclayₘₐₓMac_5::Vector{Float64}
      τ₃ₐMac_5::Vector{Float64}
      τclayΔθsrMac_5::Vector{Float64}

      τ₁ₐ_Min::Vector{Float64}
      τclay₀_Min::Vector{Float64}
      τ₂ₐ_Min::Vector{Float64}
      τclayₘₐₓ_Min::Vector{Float64}
      τ₃ₐ_Min::Vector{Float64}
      τclayΔθsr_Min::Vector{Float64}
      τ₁ₐMac_Min::Vector{Float64}
      τclay₀Mac_Min::Vector{Float64}
      τ₂ₐMac_Min::Vector{Float64}
      τclayₘₐₓMac_Min::Vector{Float64}
      τ₃ₐMac_Min::Vector{Float64}
      τclayΔθsrMac_Min::Vector{Float64}
      τ₁ₐ_1_Min::Vector{Float64}
      τclay₀_1_Min::Vector{Float64}
      τ₂ₐ_1_Min::Vector{Float64}
      τclayₘₐₓ_1_Min::Vector{Float64}
      τ₃ₐ_1_Min::Vector{Float64}
      τclayΔθsr_1_Min::Vector{Float64}

      τ₁ₐMac_1_Min::Vector{Float64}
      τclay₀Mac_1_Min::Vector{Float64}
      τ₂ₐMac_1_Min::Vector{Float64}
      τclayₘₐₓMac_1_Min::Vector{Float64}
      τ₃ₐMac_1_Min::Vector{Float64}
      τclayΔθsrMac_1_Min::Vector{Float64}

      τ₁ₐ_2_Min::Vector{Float64}
      τclay₀_2_Min::Vector{Float64}
      τ₂ₐ_2_Min::Vector{Float64}
      τclayₘₐₓ_2_Min::Vector{Float64}
      τ₃ₐ_2_Min::Vector{Float64}
      τclayΔθsr_2_Min::Vector{Float64}
      τ₁ₐMac_2_Min::Vector{Float64}
      τclay₀Mac_2_Min::Vector{Float64}
      τ₂ₐMac_2_Min::Vector{Float64}
      τclayₘₐₓMac_2_Min::Vector{Float64}
      τ₃ₐMac_2_Min::Vector{Float64}
      τclayΔθsrMac_2_Min::Vector{Float64}

      τ₁ₐ_3_Min::Vector{Float64}
      τclay₀_3_Min::Vector{Float64}
      τ₂ₐ_3_Min::Vector{Float64}
      τclayₘₐₓ_3_Min::Vector{Float64}
      τ₃ₐ_3_Min::Vector{Float64}
      τclayΔθsr_3_Min::Vector{Float64}
      τ₁ₐMac_3_Min::Vector{Float64}
      τclay₀Mac_3_Min::Vector{Float64}
      τ₂ₐMac_3_Min::Vector{Float64}
      τclayₘₐₓMac_3_Min::Vector{Float64}
      τ₃ₐMac_3_Min::Vector{Float64}
      τclayΔθsrMac_3_Min::Vector{Float64}

      τ₁ₐ_4_Min::Vector{Float64}
      τclay₀_4_Min::Vector{Float64}
      τ₂ₐ_4_Min::Vector{Float64}
      τclayₘₐₓ_4_Min::Vector{Float64}
      τ₃ₐ_4_Min::Vector{Float64}
      τclayΔθsr_4_Min::Vector{Float64}
      τ₁ₐMac_4_Min::Vector{Float64}
      τclay₀Mac_4_Min::Vector{Float64}
      τ₂ₐMac_4_Min::Vector{Float64}
      τclayₘₐₓMac_4_Min::Vector{Float64}
      τ₃ₐMac_4_Min::Vector{Float64}
      τclayΔθsrMac_4_Min::Vector{Float64}

      τ₁ₐ_5_Min::Vector{Float64}
      τclay₀_5_Min::Vector{Float64}
      τ₂ₐ_5_Min::Vector{Float64}
      τclayₘₐₓ_5_Min::Vector{Float64}
      τ₃ₐ_5_Min::Vector{Float64}
      τclayΔθsr_5_Min::Vector{Float64}
      τ₁ₐMac_5_Min::Vector{Float64}
      τclay₀Mac_5_Min::Vector{Float64}
      τ₂ₐMac_5_Min::Vector{Float64}
      τclayₘₐₓMac_5_Min::Vector{Float64}
      τ₃ₐMac_5_Min::Vector{Float64}
      τclayΔθsrMac_5_Min::Vector{Float64}

      τ₁ₐ_Max::Vector{Float64}
      τclay₀_Max::Vector{Float64}
      τ₂ₐ_Max::Vector{Float64}
      τclayₘₐₓ_Max::Vector{Float64}
      τ₃ₐ_Max::Vector{Float64}
      τclayΔθsr_Max::Vector{Float64}
      τ₁ₐMac_Max::Vector{Float64}
      τclay₀Mac_Max::Vector{Float64}
      τ₂ₐMac_Max::Vector{Float64}
      τclayₘₐₓMac_Max::Vector{Float64}
      τ₃ₐMac_Max::Vector{Float64}
      τclayΔθsrMac_Max::Vector{Float64}

      τ₁ₐ_1_Max::Vector{Float64}
      τclay₀_1_Max::Vector{Float64}
      τ₂ₐ_1_Max::Vector{Float64}
      τclayₘₐₓ_1_Max::Vector{Float64}
      τ₃ₐ_1_Max::Vector{Float64}
      τclayΔθsr_1_Max::Vector{Float64}
      τ₁ₐMac_1_Max::Vector{Float64}
      τclay₀Mac_1_Max::Vector{Float64}
      τ₂ₐMac_1_Max::Vector{Float64}
      τclayₘₐₓMac_1_Max::Vector{Float64}
      τ₃ₐMac_1_Max::Vector{Float64}
      τclayΔθsrMac_1_Max::Vector{Float64}
      τ₁ₐ_2_Max::Vector{Float64}
      τclay₀_2_Max::Vector{Float64}
      τ₂ₐ_2_Max::Vector{Float64}
      τclayₘₐₓ_2_Max::Vector{Float64}
      τ₃ₐ_2_Max::Vector{Float64}
      τclayΔθsr_2_Max::Vector{Float64}
      τ₁ₐMac_2_Max::Vector{Float64}
      τclay₀Mac_2_Max::Vector{Float64}
      τ₂ₐMac_2_Max::Vector{Float64}
      τclayₘₐₓMac_2_Max::Vector{Float64}
      τ₃ₐMac_2_Max::Vector{Float64}
      τclayΔθsrMac_2_Max::Vector{Float64}
      τ₁ₐ_3_Max::Vector{Float64}
      τclay₀_3_Max::Vector{Float64}
      τ₂ₐ_3_Max::Vector{Float64}
      τclayₘₐₓ_3_Max::Vector{Float64}
      τ₃ₐ_3_Max::Vector{Float64}
      τclayΔθsr_3_Max::Vector{Float64}
      τ₁ₐMac_3_Max::Vector{Float64}
      τclay₀Mac_3_Max::Vector{Float64}
      τ₂ₐMac_3_Max::Vector{Float64}
      τclayₘₐₓMac_3_Max::Vector{Float64}
      τ₃ₐMac_3_Max::Vector{Float64}
      τclayΔθsrMac_3_Max::Vector{Float64}
      τ₁ₐ_4_Max::Vector{Float64}
      τclay₀_4_Max::Vector{Float64}
      τ₂ₐ_4_Max::Vector{Float64}
      τclayₘₐₓ_4_Max::Vector{Float64}
      τ₃ₐ_4_Max::Vector{Float64}
      τclayΔθsr_4_Max::Vector{Float64}
      τ₁ₐMac_4_Max::Vector{Float64}
      τclay₀Mac_4_Max::Vector{Float64}
      τ₂ₐMac_4_Max::Vector{Float64}
      τclayₘₐₓMac_4_Max::Vector{Float64}
      τ₃ₐMac_4_Max::Vector{Float64}
      τclayΔθsrMac_4_Max::Vector{Float64}
      τ₁ₐ_5_Max::Vector{Float64}
      τclay₀_5_Max::Vector{Float64}
      τ₂ₐ_5_Max::Vector{Float64}
      τclayₘₐₓ_5_Max::Vector{Float64}
      τ₃ₐ_5_Max::Vector{Float64}
      τclayΔθsr_5_Max::Vector{Float64}
      τ₁ₐMac_5_Max::Vector{Float64}
      τclay₀Mac_5_Max::Vector{Float64}
      τ₂ₐMac_5_Max::Vector{Float64}
      τclayₘₐₓMac_5_Max::Vector{Float64}
      τ₃ₐMac_5_Max::Vector{Float64}
      τclayΔθsrMac_5_Max::Vector{Float64}
      
      Nse_τ::Vector{Float64}
      Rmse_τ::Vector{Float64}
      Wilmot_τ::Vector{Float64}
      Ccc_τ::Vector{Float64}
   end # mutable struct KSMODEL


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : STRUCT_KSMODEL
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function STRUCT_KSMODEL(; Nτ_Layer=Nτ_Layer::Int64)

         τ₁ₐ= fill(0.0::Float64, Nτ_Layer)
         τclay₀= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_1= fill(0.0::Float64, Nτ_Layer)
         τclay₀_1= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_1= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_1= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_1= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_1= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_1= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_1= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_1= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_1= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_1= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_1= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_2= fill(0.0::Float64, Nτ_Layer)
         τclay₀_2= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_2= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_2= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_2= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_2= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_2= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_2= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_2= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_2= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_2= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_2= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_3= fill(0.0::Float64, Nτ_Layer)
         τclay₀_3= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_3= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_3= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_3= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_3= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_3= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_3= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_3= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_3= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_3= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_3= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_4= fill(0.0::Float64, Nτ_Layer)
         τclay₀_4= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_4= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_4= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_4= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_4= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_4= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_4= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_4= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_4= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_4= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_4= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_5= fill(0.0::Float64, Nτ_Layer)
         τclay₀_5= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_5= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_5= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_5= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_5= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_5= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_5= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_5= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_5= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_5= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_5= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_1_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_1_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_1_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_2_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_2_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_2_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_3_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_3_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_3_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_4_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_4_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_4_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_5_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_5_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_5_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_1_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_1_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_1_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_2_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_2_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_2_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_3_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_3_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_3_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_4_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_4_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_4_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_5_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_5_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓ_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_5_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsr_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_5_Max= fill(0.0::Float64, Nτ_Layer)
         τclay₀Mac_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_5_Max= fill(0.0::Float64, Nτ_Layer)
         τclayₘₐₓMac_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_5_Max= fill(0.0::Float64, Nτ_Layer)
         τclayΔθsrMac_5_Max= fill(0.0::Float64, Nτ_Layer)
         

         Nse_τ      = fill(0.0::Float64, Nτ_Layer)
         Rmse_τ     = fill(0.0::Float64, Nτ_Layer)
         Wilmot_τ   = fill(0.0::Float64, Nτ_Layer)
         Ccc_τ      = fill(0.0::Float64, Nτ_Layer)

      return ksmodelτ = KSMODELτ(τ₁ₐ,τclay₀,τ₂ₐ,τclayₘₐₓ,τ₃ₐ,τclayΔθsr,τ₁ₐMac,τclay₀Mac,τ₂ₐMac,τclayₘₐₓMac,τ₃ₐMac,τclayΔθsrMac,τ₁ₐ_1,τclay₀_1,τ₂ₐ_1,τclayₘₐₓ_1,τ₃ₐ_1,τclayΔθsr_1,τ₁ₐMac_1,τclay₀Mac_1,τ₂ₐMac_1,τclayₘₐₓMac_1,τ₃ₐMac_1,τclayΔθsrMac_1,τ₁ₐ_2,τclay₀_2,τ₂ₐ_2,τclayₘₐₓ_2,τ₃ₐ_2,τclayΔθsr_2,τ₁ₐMac_2,τclay₀Mac_2,τ₂ₐMac_2,τclayₘₐₓMac_2,τ₃ₐMac_2,τclayΔθsrMac_2,τ₁ₐ_3,τclay₀_3,τ₂ₐ_3,τclayₘₐₓ_3,τ₃ₐ_3,τclayΔθsr_3,τ₁ₐMac_3,τclay₀Mac_3,τ₂ₐMac_3,τclayₘₐₓMac_3,τ₃ₐMac_3,τclayΔθsrMac_3,τ₁ₐ_4,τclay₀_4,τ₂ₐ_4,τclayₘₐₓ_4,τ₃ₐ_4,τclayΔθsr_4,τ₁ₐMac_4,τclay₀Mac_4,τ₂ₐMac_4,τclayₘₐₓMac_4,τ₃ₐMac_4,τclayΔθsrMac_4,τ₁ₐ_5,τclay₀_5,τ₂ₐ_5,τclayₘₐₓ_5,τ₃ₐ_5,τclayΔθsr_5,τ₁ₐMac_5,τclay₀Mac_5,τ₂ₐMac_5,τclayₘₐₓMac_5,τ₃ₐMac_5,τclayΔθsrMac_5,τ₁ₐ_Min,τclay₀_Min,τ₂ₐ_Min,τclayₘₐₓ_Min,τ₃ₐ_Min,τclayΔθsr_Min,τ₁ₐMac_Min,τclay₀Mac_Min,τ₂ₐMac_Min,τclayₘₐₓMac_Min,τ₃ₐMac_Min,τclayΔθsrMac_Min,τ₁ₐ_1_Min,τclay₀_1_Min,τ₂ₐ_1_Min,τclayₘₐₓ_1_Min,τ₃ₐ_1_Min,τclayΔθsr_1_Min,τ₁ₐMac_1_Min,τclay₀Mac_1_Min,τ₂ₐMac_1_Min,τclayₘₐₓMac_1_Min,τ₃ₐMac_1_Min,τclayΔθsrMac_1_Min,τ₁ₐ_2_Min,τclay₀_2_Min,τ₂ₐ_2_Min,τclayₘₐₓ_2_Min,τ₃ₐ_2_Min,τclayΔθsr_2_Min,τ₁ₐMac_2_Min,τclay₀Mac_2_Min,τ₂ₐMac_2_Min,τclayₘₐₓMac_2_Min,τ₃ₐMac_2_Min,τclayΔθsrMac_2_Min,τ₁ₐ_3_Min,τclay₀_3_Min,τ₂ₐ_3_Min,τclayₘₐₓ_3_Min,τ₃ₐ_3_Min,τclayΔθsr_3_Min,τ₁ₐMac_3_Min,τclay₀Mac_3_Min,τ₂ₐMac_3_Min,τclayₘₐₓMac_3_Min,τ₃ₐMac_3_Min,τclayΔθsrMac_3_Min,τ₁ₐ_4_Min,τclay₀_4_Min,τ₂ₐ_4_Min,τclayₘₐₓ_4_Min,τ₃ₐ_4_Min,τclayΔθsr_4_Min,τ₁ₐMac_4_Min,τclay₀Mac_4_Min,τ₂ₐMac_4_Min,τclayₘₐₓMac_4_Min,τ₃ₐMac_4_Min,τclayΔθsrMac_4_Min,τ₁ₐ_5_Min,τclay₀_5_Min,τ₂ₐ_5_Min,τclayₘₐₓ_5_Min,τ₃ₐ_5_Min,τclayΔθsr_5_Min,τ₁ₐMac_5_Min,τclay₀Mac_5_Min,τ₂ₐMac_5_Min,τclayₘₐₓMac_5_Min,τ₃ₐMac_5_Min,τclayΔθsrMac_5_Min,τ₁ₐ_Max,τclay₀_Max,τ₂ₐ_Max,τclayₘₐₓ_Max,τ₃ₐ_Max,τclayΔθsr_Max,τ₁ₐMac_Max,τclay₀Mac_Max,τ₂ₐMac_Max,τclayₘₐₓMac_Max,τ₃ₐMac_Max,τclayΔθsrMac_Max,τ₁ₐ_1_Max,τclay₀_1_Max,τ₂ₐ_1_Max,τclayₘₐₓ_1_Max,τ₃ₐ_1_Max,τclayΔθsr_1_Max,τ₁ₐMac_1_Max,τclay₀Mac_1_Max,τ₂ₐMac_1_Max,τclayₘₐₓMac_1_Max,τ₃ₐMac_1_Max,τclayΔθsrMac_1_Max,τ₁ₐ_2_Max,τclay₀_2_Max,τ₂ₐ_2_Max,τclayₘₐₓ_2_Max,τ₃ₐ_2_Max,τclayΔθsr_2_Max,τ₁ₐMac_2_Max,τclay₀Mac_2_Max,τ₂ₐMac_2_Max,τclayₘₐₓMac_2_Max,τ₃ₐMac_2_Max,τclayΔθsrMac_2_Max,τ₁ₐ_3_Max,τclay₀_3_Max,τ₂ₐ_3_Max,τclayₘₐₓ_3_Max,τ₃ₐ_3_Max,τclayΔθsr_3_Max,τ₁ₐMac_3_Max,τclay₀Mac_3_Max,τ₂ₐMac_3_Max,τclayₘₐₓMac_3_Max,τ₃ₐMac_3_Max,τclayΔθsrMac_3_Max,τ₁ₐ_4_Max,τclay₀_4_Max,τ₂ₐ_4_Max,τclayₘₐₓ_4_Max,τ₃ₐ_4_Max,τclayΔθsr_4_Max,τ₁ₐMac_4_Max,τclay₀Mac_4_Max,τ₂ₐMac_4_Max,τclayₘₐₓMac_4_Max,τ₃ₐMac_4_Max,τclayΔθsrMac_4_Max,τ₁ₐ_5_Max,τclay₀_5_Max,τ₂ₐ_5_Max,τclayₘₐₓ_5_Max,τ₃ₐ_5_Max,τclayΔθsr_5_Max,τ₁ₐMac_5_Max,τclay₀Mac_5_Max,τ₂ₐMac_5_Max,τclayₘₐₓMac_5_Max,τ₃ₐMac_5_Max,τclayΔθsrMac_5_Max, Nse_τ,Rmse_τ,Wilmot_τ,Ccc_τ)

      end  # function: STRUCT_KSMODEL

   end  # module: ksModel
# ............................................................
