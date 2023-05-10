# =============================================================
#		module: hydroSmooth
# =============================================================
module hydroSmooth

   import ..tool, ..hydroStruct
   export HYDROHORIZON_2_HYDRO_SMOOTENING, DISCRETISATION_SMOOTENING!

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : SMOOTENING_DISCRETISATION
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function DISCRETISATION_SMOOTENING!(🎏_θΨini::Symbol, iScenario::Int64, N_SoilLayer::Int64, optionHypix, paramHypix, pathInputHypix, Zlayer::Vector{Float64}, θini_or_Ψini::Vector{Float64})
            
         ΔZfine_⬓ = paramHypix.ΔZfine * 0.5

         # Read hydro parameters
            hydro_Smooth    = hydroStruct.HYDROSTRUCT(optionHypix, N_SoilLayer)
            hydro_Smooth, ~ = tool.readWrite.READ_STRUCT_SIMPLE(hydro_Smooth, pathInputHypix.HydroInput[iScenario])

         # New layering
            Layer₀        = Float64[]
            θini_or_Ψini₀ = Float64[]
            Zlayer₀       = Float64[]

         for iZ=1:N_SoilLayer-1
            append!(Layer₀, iZ)
            append!(θini_or_Ψini₀, θini_or_Ψini[iZ])
            append!(Zlayer₀, Zlayer[iZ])

            # Determening if smootening layer is going to put
            if (abs((hydro_Smooth.θs[iZ] - hydro_Smooth.θr[iZ]) - (hydro_Smooth.θs[iZ+1] - hydro_Smooth.θr[iZ+1])) ≥  paramHypix.ΔθSθr_Smooth) || ((abs(log10(hydro_Smooth.Ks[iZ]) - log10(hydro_Smooth.Ks[iZ+1]))) ≥ paramHypix.ΔKsLog_Smooth) && (Zlayer[iZ] + ΔZfine_⬓ < Zlayer[iZ+1] - ΔZfine_⬓)

               printstyled("    === Smootening layer:  " , iZ, "   === \n"; color=:green) 

               # Adding smootening Layer
                  append!(Layer₀, iZ + 0.5)

               # Adding smootening θini_or_Ψini₀
                  if 🎏_θΨini == :θini
                     append!(θini_or_Ψini₀, (θini_or_Ψini[iZ] + θini_or_Ψini[iZ+1]) * 0.5)
                  else
                     append!(θini_or_Ψini₀, exp((log(θini_or_Ψini[iZ]) + log(θini_or_Ψini[iZ+1])) * 0.5))
                  end

               # Adding Zlayer
                  Zlayer₀[end] =  Zlayer₀[end] - ΔZfine_⬓
                  append!(Zlayer₀, Zlayer[iZ] + ΔZfine_⬓) 
            end #if
         end #  for iZ=1:N_SoilLayer-1
         append!(θini_or_Ψini₀, θini_or_Ψini[N_SoilLayer])
         append!(Layer₀, N_SoilLayer)
         append!(Zlayer₀, Zlayer[N_SoilLayer])
         N_Layer = length(Layer₀)

      return Layer₀, N_Layer, Zlayer₀, θini_or_Ψini₀  
      end # function SMOOTENING_DISCRETISATION
   #-------------------------------------------------------------------------


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : HYDROHORIZON_2_HYDRO_SMOOTENING
   #                 ONLY FOR KOSUGI MODEL
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function HYDROHORIZON_2_HYDRO_SMOOTENING(discret, hydroHorizon, Layer_Smooth, optionHypix)

         Kₛ_AVER(Kₛ₁, Kₛ₂, w) = exp(w *log(Kₛ₁) + (1.0 - w) *log(Kₛ₂))
      
         σ_AVER(σ₁, σ₂, w) = √(w * σ₁ ^ 2.0 + (1.0 - w) * σ₂ ^ 2.0)

         Ψₘ_AVER(Ψₘ₁, Ψₘ₂, w) =exp((w *log(Ψₘ₁) + (1.0 - w) * log(Ψₘ₂)))

         θₛ_AVER(θₛ₁, θₛ₂, w) = (w * θₛ₁ + (1.0 - w) * θₛ₂)

         θᵣ_AVER(θᵣ₁, θᵣ₂, w) = (w * θᵣ₁ + (1.0 - w) * θᵣ₂)

         θₛMacMat_AVER(θₛMacMat₁, θₛMacMat₂, w) = (w * θₛMacMat₁ + (1.0 - w) * θₛMacMat₂)

         N_Layer = length(Layer_Smooth)
         
         # Hydraulic parameters with smootening layers
            hydro_Smooth = hydroStruct.HYDROSTRUCT(optionHypix, N_Layer)

         # INITIALISING THE STRUCTURE
         # looping through every fieldnames of the structure
            FieldName_Array = propertynames(hydroHorizon)

            Vector_FieldName = fill(0.0::Float64, N_Layer, length(FieldName_Array))

         iZ = 1
         for iSmooth = 1:N_Layer
            # Determening the value of the other hydraulic parameters
               iFieldName = 1
               for FieldName in FieldName_Array
                  if  occursin("_Min", String(FieldName)) || occursin("_Max", String(FieldName))
                     Vector_FieldName[iSmooth, iFieldName] = getfield(hydroHorizon, FieldName)[iZ]	
                     setfield!(hydro_Smooth, Symbol(FieldName), Vector_FieldName[:, iFieldName])
                  end
                  iFieldName += 1
               end #  for FieldName in FieldName_Array

            # Determening if it is a smootening layer is detected as it number is i + 0.5
            if Layer_Smooth[iSmooth] - floor(Layer_Smooth[iSmooth]) > 0.4

               iZ = Int64(Layer_Smooth[iSmooth+1])

               hydro_Smooth.Ks[iSmooth]       = Kₛ_AVER(hydroHorizon.Ks[iZ], hydroHorizon.Ks[iZ-1], discret.ΔZ_W[iSmooth])

               hydro_Smooth.ΨmMac[iSmooth]    = Ψₘ_AVER(hydroHorizon.ΨmMac[iZ], hydroHorizon.ΨmMac[iZ-1], 0.5)

               hydro_Smooth.Ψm[iSmooth]       = Ψₘ_AVER(hydroHorizon.Ψm[iZ], hydroHorizon.Ψm[iZ-1], 0.5)

               hydro_Smooth.θr[iSmooth]       = θᵣ_AVER(hydroHorizon.θr[iZ], hydroHorizon.θr[iZ-1], 0.5)

               hydro_Smooth.θsMacMat[iSmooth] = θₛMacMat_AVER(hydroHorizon.θsMacMat[iZ], hydroHorizon.θsMacMat[iZ-1], 0.5)

               hydro_Smooth.θs[iSmooth]       = θₛ_AVER(hydroHorizon.θs[iZ], hydroHorizon.θs[iZ-1], 0.5)

               hydro_Smooth.σMac[iSmooth]     = σ_AVER(hydroHorizon.σMac[iZ], hydroHorizon.σMac[iZ-1], 0.5)

               hydro_Smooth.σ[iSmooth]        = σ_AVER(hydroHorizon.σ[iZ], hydroHorizon.σ[iZ-1], 0.5)

               hydro_Smooth.θsMacMat_ƞ[iSmooth] = min((hydroHorizon.θsMacMat[iZ] - hydroHorizon.θr[iZ]) / (hydroHorizon.θs[iZ] - hydroHorizon.θr[iZ]), 1.0)

            else
               iZ = Int64(Layer_Smooth[iSmooth])

               hydro_Smooth.Ks[iSmooth]       = hydroHorizon.Ks[iZ]

               hydro_Smooth.ΨmMac[iSmooth]    = hydroHorizon.ΨmMac[iZ]

               hydro_Smooth.Ψm[iSmooth]       = hydroHorizon.Ψm[iZ]

               hydro_Smooth.θr[iSmooth]       = hydroHorizon.θr[iZ]

               hydro_Smooth.θsMacMat[iSmooth] = hydroHorizon.θsMacMat[iZ]

               hydro_Smooth.θs[iSmooth]       = hydroHorizon.θs[iZ]

               hydro_Smooth.σMac[iSmooth]     = hydroHorizon.σMac[iZ]

               hydro_Smooth.σ[iSmooth]        = hydroHorizon.σ[iZ]

               hydro_Smooth.θsMacMat_ƞ[iSmooth] = min((hydroHorizon.θsMacMat[iZ] - hydroHorizon.θr[iZ]) / (hydroHorizon.θs[iZ] - hydroHorizon.θr[iZ]), 1.0)
              
            end #  if Layer_Smooth[iSmooth] - floor(Layer_Smooth[iSmooth]) > 0.4

         end #  for iSmooth = 1:N_Layer

      return hydro_Smooth  
      end # function HYDROHORIZON_2_HYDRO_SMOOTENING
   #-------------------------------------------------------------------------

   
end  # module: hydroSmooth
# ............................................................