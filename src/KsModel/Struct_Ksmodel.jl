# =============================================================
#		module: ksModel
# =============================================================
module ksModel

   Base.@kwdef mutable struct KSMODELτ
      τ₁ₐ::Vector{Float64}
      τclay::Vector{Float64}
      τ₂ₐ::Vector{Float64}
      τ₂ᵦ::Vector{Float64}
      τ₃ₐ::Vector{Float64}
      τ₃ᵦ::Vector{Float64}
      τ₁ₐMac::Vector{Float64}
      τclayMac::Vector{Float64}
      τ₂ₐMac::Vector{Float64}
      τ₂ᵦMac::Vector{Float64}
      τ₃ₐMac::Vector{Float64}
      τ₃ᵦMac::Vector{Float64}

      τ₁ₐ_1::Vector{Float64}
      τclay_1::Vector{Float64}
      τ₂ₐ_1::Vector{Float64}
      τ₂ᵦ_1::Vector{Float64}
      τ₃ₐ_1::Vector{Float64}
      τ₃ᵦ_1::Vector{Float64}
      τ₁ₐMac_1::Vector{Float64}
      τclayMac_1::Vector{Float64}
      τ₂ₐMac_1::Vector{Float64}
      τ₂ᵦMac_1::Vector{Float64}
      τ₃ₐMac_1::Vector{Float64}
      τ₃ᵦMac_1::Vector{Float64}

      τ₁ₐ_2::Vector{Float64}
      τclay_2::Vector{Float64}
      τ₂ₐ_2::Vector{Float64}
      τ₂ᵦ_2::Vector{Float64}
      τ₃ₐ_2::Vector{Float64}
      τ₃ᵦ_2::Vector{Float64}
      τ₁ₐMac_2::Vector{Float64}
      τclayMac_2::Vector{Float64}
      τ₂ₐMac_2::Vector{Float64}
      τ₂ᵦMac_2::Vector{Float64}
      τ₃ₐMac_2::Vector{Float64}
      τ₃ᵦMac_2::Vector{Float64}

      τ₁ₐ_3::Vector{Float64}
      τclay_3::Vector{Float64}
      τ₂ₐ_3::Vector{Float64}
      τ₂ᵦ_3::Vector{Float64}
      τ₃ₐ_3::Vector{Float64}
      τ₃ᵦ_3::Vector{Float64}
      τ₁ₐMac_3::Vector{Float64}
      τclayMac_3::Vector{Float64}
      τ₂ₐMac_3::Vector{Float64}
      τ₂ᵦMac_3::Vector{Float64}
      τ₃ₐMac_3::Vector{Float64}
      τ₃ᵦMac_3::Vector{Float64}

      τ₁ₐ_4::Vector{Float64}
      τclay_4::Vector{Float64}
      τ₂ₐ_4::Vector{Float64}
      τ₂ᵦ_4::Vector{Float64}
      τ₃ₐ_4::Vector{Float64}
      τ₃ᵦ_4::Vector{Float64}
      τ₁ₐMac_4::Vector{Float64}
      τclayMac_4::Vector{Float64}
      τ₂ₐMac_4::Vector{Float64}
      τ₂ᵦMac_4::Vector{Float64}
      τ₃ₐMac_4::Vector{Float64}
      τ₃ᵦMac_4::Vector{Float64}

      τ₁ₐ_5::Vector{Float64}
      τclay_5::Vector{Float64}
      τ₂ₐ_5::Vector{Float64}
      τ₂ᵦ_5::Vector{Float64}
      τ₃ₐ_5::Vector{Float64}
      τ₃ᵦ_5::Vector{Float64}
      τ₁ₐMac_5::Vector{Float64}
      τclayMac_5::Vector{Float64}
      τ₂ₐMac_5::Vector{Float64}
      τ₂ᵦMac_5::Vector{Float64}
      τ₃ₐMac_5::Vector{Float64}
      τ₃ᵦMac_5::Vector{Float64}

      τ₁ₐ_Min::Vector{Float64}
      τclay_Min::Vector{Float64}
      τ₂ₐ_Min::Vector{Float64}
      τ₂ᵦ_Min::Vector{Float64}
      τ₃ₐ_Min::Vector{Float64}
      τ₃ᵦ_Min::Vector{Float64}
      τ₁ₐMac_Min::Vector{Float64}
      τclayMac_Min::Vector{Float64}
      τ₂ₐMac_Min::Vector{Float64}
      τ₂ᵦMac_Min::Vector{Float64}
      τ₃ₐMac_Min::Vector{Float64}
      τ₃ᵦMac_Min::Vector{Float64}
      τ₁ₐ_1_Min::Vector{Float64}
      τclay_1_Min::Vector{Float64}
      τ₂ₐ_1_Min::Vector{Float64}
      τ₂ᵦ_1_Min::Vector{Float64}
      τ₃ₐ_1_Min::Vector{Float64}
      τ₃ᵦ_1_Min::Vector{Float64}

      τ₁ₐMac_1_Min::Vector{Float64}
      τclayMac_1_Min::Vector{Float64}
      τ₂ₐMac_1_Min::Vector{Float64}
      τ₂ᵦMac_1_Min::Vector{Float64}
      τ₃ₐMac_1_Min::Vector{Float64}
      τ₃ᵦMac_1_Min::Vector{Float64}

      τ₁ₐ_2_Min::Vector{Float64}
      τclay_2_Min::Vector{Float64}
      τ₂ₐ_2_Min::Vector{Float64}
      τ₂ᵦ_2_Min::Vector{Float64}
      τ₃ₐ_2_Min::Vector{Float64}
      τ₃ᵦ_2_Min::Vector{Float64}
      τ₁ₐMac_2_Min::Vector{Float64}
      τclayMac_2_Min::Vector{Float64}
      τ₂ₐMac_2_Min::Vector{Float64}
      τ₂ᵦMac_2_Min::Vector{Float64}
      τ₃ₐMac_2_Min::Vector{Float64}
      τ₃ᵦMac_2_Min::Vector{Float64}

      τ₁ₐ_3_Min::Vector{Float64}
      τclay_3_Min::Vector{Float64}
      τ₂ₐ_3_Min::Vector{Float64}
      τ₂ᵦ_3_Min::Vector{Float64}
      τ₃ₐ_3_Min::Vector{Float64}
      τ₃ᵦ_3_Min::Vector{Float64}
      τ₁ₐMac_3_Min::Vector{Float64}
      τclayMac_3_Min::Vector{Float64}
      τ₂ₐMac_3_Min::Vector{Float64}
      τ₂ᵦMac_3_Min::Vector{Float64}
      τ₃ₐMac_3_Min::Vector{Float64}
      τ₃ᵦMac_3_Min::Vector{Float64}

      τ₁ₐ_4_Min::Vector{Float64}
      τclay_4_Min::Vector{Float64}
      τ₂ₐ_4_Min::Vector{Float64}
      τ₂ᵦ_4_Min::Vector{Float64}
      τ₃ₐ_4_Min::Vector{Float64}
      τ₃ᵦ_4_Min::Vector{Float64}
      τ₁ₐMac_4_Min::Vector{Float64}
      τclayMac_4_Min::Vector{Float64}
      τ₂ₐMac_4_Min::Vector{Float64}
      τ₂ᵦMac_4_Min::Vector{Float64}
      τ₃ₐMac_4_Min::Vector{Float64}
      τ₃ᵦMac_4_Min::Vector{Float64}

      τ₁ₐ_5_Min::Vector{Float64}
      τclay_5_Min::Vector{Float64}
      τ₂ₐ_5_Min::Vector{Float64}
      τ₂ᵦ_5_Min::Vector{Float64}
      τ₃ₐ_5_Min::Vector{Float64}
      τ₃ᵦ_5_Min::Vector{Float64}
      τ₁ₐMac_5_Min::Vector{Float64}
      τclayMac_5_Min::Vector{Float64}
      τ₂ₐMac_5_Min::Vector{Float64}
      τ₂ᵦMac_5_Min::Vector{Float64}
      τ₃ₐMac_5_Min::Vector{Float64}
      τ₃ᵦMac_5_Min::Vector{Float64}

      τ₁ₐ_Max::Vector{Float64}
      τclay_Max::Vector{Float64}
      τ₂ₐ_Max::Vector{Float64}
      τ₂ᵦ_Max::Vector{Float64}
      τ₃ₐ_Max::Vector{Float64}
      τ₃ᵦ_Max::Vector{Float64}
      τ₁ₐMac_Max::Vector{Float64}
      τclayMac_Max::Vector{Float64}
      τ₂ₐMac_Max::Vector{Float64}
      τ₂ᵦMac_Max::Vector{Float64}
      τ₃ₐMac_Max::Vector{Float64}
      τ₃ᵦMac_Max::Vector{Float64}

      τ₁ₐ_1_Max::Vector{Float64}
      τclay_1_Max::Vector{Float64}
      τ₂ₐ_1_Max::Vector{Float64}
      τ₂ᵦ_1_Max::Vector{Float64}
      τ₃ₐ_1_Max::Vector{Float64}
      τ₃ᵦ_1_Max::Vector{Float64}
      τ₁ₐMac_1_Max::Vector{Float64}
      τclayMac_1_Max::Vector{Float64}
      τ₂ₐMac_1_Max::Vector{Float64}
      τ₂ᵦMac_1_Max::Vector{Float64}
      τ₃ₐMac_1_Max::Vector{Float64}
      τ₃ᵦMac_1_Max::Vector{Float64}
      τ₁ₐ_2_Max::Vector{Float64}
      τclay_2_Max::Vector{Float64}
      τ₂ₐ_2_Max::Vector{Float64}
      τ₂ᵦ_2_Max::Vector{Float64}
      τ₃ₐ_2_Max::Vector{Float64}
      τ₃ᵦ_2_Max::Vector{Float64}
      τ₁ₐMac_2_Max::Vector{Float64}
      τclayMac_2_Max::Vector{Float64}
      τ₂ₐMac_2_Max::Vector{Float64}
      τ₂ᵦMac_2_Max::Vector{Float64}
      τ₃ₐMac_2_Max::Vector{Float64}
      τ₃ᵦMac_2_Max::Vector{Float64}
      τ₁ₐ_3_Max::Vector{Float64}
      τclay_3_Max::Vector{Float64}
      τ₂ₐ_3_Max::Vector{Float64}
      τ₂ᵦ_3_Max::Vector{Float64}
      τ₃ₐ_3_Max::Vector{Float64}
      τ₃ᵦ_3_Max::Vector{Float64}
      τ₁ₐMac_3_Max::Vector{Float64}
      τclayMac_3_Max::Vector{Float64}
      τ₂ₐMac_3_Max::Vector{Float64}
      τ₂ᵦMac_3_Max::Vector{Float64}
      τ₃ₐMac_3_Max::Vector{Float64}
      τ₃ᵦMac_3_Max::Vector{Float64}
      τ₁ₐ_4_Max::Vector{Float64}
      τclay_4_Max::Vector{Float64}
      τ₂ₐ_4_Max::Vector{Float64}
      τ₂ᵦ_4_Max::Vector{Float64}
      τ₃ₐ_4_Max::Vector{Float64}
      τ₃ᵦ_4_Max::Vector{Float64}
      τ₁ₐMac_4_Max::Vector{Float64}
      τclayMac_4_Max::Vector{Float64}
      τ₂ₐMac_4_Max::Vector{Float64}
      τ₂ᵦMac_4_Max::Vector{Float64}
      τ₃ₐMac_4_Max::Vector{Float64}
      τ₃ᵦMac_4_Max::Vector{Float64}
      τ₁ₐ_5_Max::Vector{Float64}
      τclay_5_Max::Vector{Float64}
      τ₂ₐ_5_Max::Vector{Float64}
      τ₂ᵦ_5_Max::Vector{Float64}
      τ₃ₐ_5_Max::Vector{Float64}
      τ₃ᵦ_5_Max::Vector{Float64}
      τ₁ₐMac_5_Max::Vector{Float64}
      τclayMac_5_Max::Vector{Float64}
      τ₂ₐMac_5_Max::Vector{Float64}
      τ₂ᵦMac_5_Max::Vector{Float64}
      τ₃ₐMac_5_Max::Vector{Float64}
      τ₃ᵦMac_5_Max::Vector{Float64}
      
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
         τclay= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac= fill(0.0::Float64, Nτ_Layer)
         τclayMac= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_1= fill(0.0::Float64, Nτ_Layer)
         τclay_1= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_1= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_1= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_1= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_1= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_1= fill(0.0::Float64, Nτ_Layer)
         τclayMac_1= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_1= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_1= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_1= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_1= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_2= fill(0.0::Float64, Nτ_Layer)
         τclay_2= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_2= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_2= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_2= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_2= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_2= fill(0.0::Float64, Nτ_Layer)
         τclayMac_2= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_2= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_2= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_2= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_2= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_3= fill(0.0::Float64, Nτ_Layer)
         τclay_3= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_3= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_3= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_3= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_3= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_3= fill(0.0::Float64, Nτ_Layer)
         τclayMac_3= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_3= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_3= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_3= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_3= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_4= fill(0.0::Float64, Nτ_Layer)
         τclay_4= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_4= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_4= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_4= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_4= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_4= fill(0.0::Float64, Nτ_Layer)
         τclayMac_4= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_4= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_4= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_4= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_4= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_5= fill(0.0::Float64, Nτ_Layer)
         τclay_5= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_5= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_5= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_5= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_5= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_5= fill(0.0::Float64, Nτ_Layer)
         τclayMac_5= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_5= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_5= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_5= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_5= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_Min= fill(0.0::Float64, Nτ_Layer)
         τclay_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_Min= fill(0.0::Float64, Nτ_Layer)
         τclayMac_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_1_Min= fill(0.0::Float64, Nτ_Layer)
         τclay_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τclayMac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_1_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_2_Min= fill(0.0::Float64, Nτ_Layer)
         τclay_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τclayMac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_2_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_3_Min= fill(0.0::Float64, Nτ_Layer)
         τclay_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τclayMac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_3_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_4_Min= fill(0.0::Float64, Nτ_Layer)
         τclay_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τclayMac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_4_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_5_Min= fill(0.0::Float64, Nτ_Layer)
         τclay_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τclayMac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_5_Min= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_Max= fill(0.0::Float64, Nτ_Layer)
         τclay_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_Max= fill(0.0::Float64, Nτ_Layer)
         τclayMac_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_1_Max= fill(0.0::Float64, Nτ_Layer)
         τclay_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τclayMac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_1_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_2_Max= fill(0.0::Float64, Nτ_Layer)
         τclay_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τclayMac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_2_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_3_Max= fill(0.0::Float64, Nτ_Layer)
         τclay_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τclayMac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_3_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_4_Max= fill(0.0::Float64, Nτ_Layer)
         τclay_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τclayMac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_4_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐ_5_Max= fill(0.0::Float64, Nτ_Layer)
         τclay_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐ_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦ_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐ_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦ_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₁ₐMac_5_Max= fill(0.0::Float64, Nτ_Layer)
         τclayMac_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ₐMac_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₂ᵦMac_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ₐMac_5_Max= fill(0.0::Float64, Nτ_Layer)
         τ₃ᵦMac_5_Max= fill(0.0::Float64, Nτ_Layer)
         

         Nse_τ      = fill(0.0::Float64, Nτ_Layer)
         Rmse_τ     = fill(0.0::Float64, Nτ_Layer)
         Wilmot_τ   = fill(0.0::Float64, Nτ_Layer)
         Ccc_τ      = fill(0.0::Float64, Nτ_Layer)

      return ksmodelτ = KSMODELτ(τ₁ₐ,τclay,τ₂ₐ,τ₂ᵦ,τ₃ₐ,τ₃ᵦ,τ₁ₐMac,τclayMac,τ₂ₐMac,τ₂ᵦMac,τ₃ₐMac,τ₃ᵦMac,τ₁ₐ_1,τclay_1,τ₂ₐ_1,τ₂ᵦ_1,τ₃ₐ_1,τ₃ᵦ_1,τ₁ₐMac_1,τclayMac_1,τ₂ₐMac_1,τ₂ᵦMac_1,τ₃ₐMac_1,τ₃ᵦMac_1,τ₁ₐ_2,τclay_2,τ₂ₐ_2,τ₂ᵦ_2,τ₃ₐ_2,τ₃ᵦ_2,τ₁ₐMac_2,τclayMac_2,τ₂ₐMac_2,τ₂ᵦMac_2,τ₃ₐMac_2,τ₃ᵦMac_2,τ₁ₐ_3,τclay_3,τ₂ₐ_3,τ₂ᵦ_3,τ₃ₐ_3,τ₃ᵦ_3,τ₁ₐMac_3,τclayMac_3,τ₂ₐMac_3,τ₂ᵦMac_3,τ₃ₐMac_3,τ₃ᵦMac_3,τ₁ₐ_4,τclay_4,τ₂ₐ_4,τ₂ᵦ_4,τ₃ₐ_4,τ₃ᵦ_4,τ₁ₐMac_4,τclayMac_4,τ₂ₐMac_4,τ₂ᵦMac_4,τ₃ₐMac_4,τ₃ᵦMac_4,τ₁ₐ_5,τclay_5,τ₂ₐ_5,τ₂ᵦ_5,τ₃ₐ_5,τ₃ᵦ_5,τ₁ₐMac_5,τclayMac_5,τ₂ₐMac_5,τ₂ᵦMac_5,τ₃ₐMac_5,τ₃ᵦMac_5,τ₁ₐ_Min,τclay_Min,τ₂ₐ_Min,τ₂ᵦ_Min,τ₃ₐ_Min,τ₃ᵦ_Min,τ₁ₐMac_Min,τclayMac_Min,τ₂ₐMac_Min,τ₂ᵦMac_Min,τ₃ₐMac_Min,τ₃ᵦMac_Min,τ₁ₐ_1_Min,τclay_1_Min,τ₂ₐ_1_Min,τ₂ᵦ_1_Min,τ₃ₐ_1_Min,τ₃ᵦ_1_Min,τ₁ₐMac_1_Min,τclayMac_1_Min,τ₂ₐMac_1_Min,τ₂ᵦMac_1_Min,τ₃ₐMac_1_Min,τ₃ᵦMac_1_Min,τ₁ₐ_2_Min,τclay_2_Min,τ₂ₐ_2_Min,τ₂ᵦ_2_Min,τ₃ₐ_2_Min,τ₃ᵦ_2_Min,τ₁ₐMac_2_Min,τclayMac_2_Min,τ₂ₐMac_2_Min,τ₂ᵦMac_2_Min,τ₃ₐMac_2_Min,τ₃ᵦMac_2_Min,τ₁ₐ_3_Min,τclay_3_Min,τ₂ₐ_3_Min,τ₂ᵦ_3_Min,τ₃ₐ_3_Min,τ₃ᵦ_3_Min,τ₁ₐMac_3_Min,τclayMac_3_Min,τ₂ₐMac_3_Min,τ₂ᵦMac_3_Min,τ₃ₐMac_3_Min,τ₃ᵦMac_3_Min,τ₁ₐ_4_Min,τclay_4_Min,τ₂ₐ_4_Min,τ₂ᵦ_4_Min,τ₃ₐ_4_Min,τ₃ᵦ_4_Min,τ₁ₐMac_4_Min,τclayMac_4_Min,τ₂ₐMac_4_Min,τ₂ᵦMac_4_Min,τ₃ₐMac_4_Min,τ₃ᵦMac_4_Min,τ₁ₐ_5_Min,τclay_5_Min,τ₂ₐ_5_Min,τ₂ᵦ_5_Min,τ₃ₐ_5_Min,τ₃ᵦ_5_Min,τ₁ₐMac_5_Min,τclayMac_5_Min,τ₂ₐMac_5_Min,τ₂ᵦMac_5_Min,τ₃ₐMac_5_Min,τ₃ᵦMac_5_Min,τ₁ₐ_Max,τclay_Max,τ₂ₐ_Max,τ₂ᵦ_Max,τ₃ₐ_Max,τ₃ᵦ_Max,τ₁ₐMac_Max,τclayMac_Max,τ₂ₐMac_Max,τ₂ᵦMac_Max,τ₃ₐMac_Max,τ₃ᵦMac_Max,τ₁ₐ_1_Max,τclay_1_Max,τ₂ₐ_1_Max,τ₂ᵦ_1_Max,τ₃ₐ_1_Max,τ₃ᵦ_1_Max,τ₁ₐMac_1_Max,τclayMac_1_Max,τ₂ₐMac_1_Max,τ₂ᵦMac_1_Max,τ₃ₐMac_1_Max,τ₃ᵦMac_1_Max,τ₁ₐ_2_Max,τclay_2_Max,τ₂ₐ_2_Max,τ₂ᵦ_2_Max,τ₃ₐ_2_Max,τ₃ᵦ_2_Max,τ₁ₐMac_2_Max,τclayMac_2_Max,τ₂ₐMac_2_Max,τ₂ᵦMac_2_Max,τ₃ₐMac_2_Max,τ₃ᵦMac_2_Max,τ₁ₐ_3_Max,τclay_3_Max,τ₂ₐ_3_Max,τ₂ᵦ_3_Max,τ₃ₐ_3_Max,τ₃ᵦ_3_Max,τ₁ₐMac_3_Max,τclayMac_3_Max,τ₂ₐMac_3_Max,τ₂ᵦMac_3_Max,τ₃ₐMac_3_Max,τ₃ᵦMac_3_Max,τ₁ₐ_4_Max,τclay_4_Max,τ₂ₐ_4_Max,τ₂ᵦ_4_Max,τ₃ₐ_4_Max,τ₃ᵦ_4_Max,τ₁ₐMac_4_Max,τclayMac_4_Max,τ₂ₐMac_4_Max,τ₂ᵦMac_4_Max,τ₃ₐMac_4_Max,τ₃ᵦMac_4_Max,τ₁ₐ_5_Max,τclay_5_Max,τ₂ₐ_5_Max,τ₂ᵦ_5_Max,τ₃ₐ_5_Max,τ₃ᵦ_5_Max,τ₁ₐMac_5_Max,τclayMac_5_Max,τ₂ₐMac_5_Max,τ₂ᵦMac_5_Max,τ₃ₐMac_5_Max,τ₃ᵦMac_5_Max, Nse_τ,Rmse_τ,Wilmot_τ,Ccc_τ)

      end  # function: STRUCT_KSMODEL

   end  # module: ksModel
# ............................................................
