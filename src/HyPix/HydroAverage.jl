# =============================================================
#		module: hydroAver
# =============================================================
module hydroAver

   import Yeppp

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : HYDRO_AVERAGE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function HYDRO_AVERAGE(hydro, Nz)

      Kₛ_AVER(Kₛ₁, Kₛ₂, w) = Yeppp.exp(w * Yeppp.log(Kₛ₁) + (1.0 - w) * Yeppp.log(Kₛ₂))
      
      σ_AVER(σ₁, σ₂, w) = √(w * σ₁ ^ 2.0 + (1.0 - w) * σ₂ ^ 2.0)

      Ψₘ_AVER(Ψₘ₁, Ψₘ₂, w) = Yeppp.exp((w * Yeppp.log(Ψₘ₁) + (1.0 - w) * Yeppp.log(Ψₘ₂)) * 0.5)

      θₛ_AVER(θₛ₁, θₛ₂, w) = (w * θₛ₁ + (1.0 - w) * θₛ₂)

      θᵣ_AVER(θᵣ₁, θᵣ₂, w) = (w * θᵣ₁ + (1.0 - w) * θᵣ₂)

      θₛMac_AVER(θₛMac₁, θₛMac₂, w) = (w * θₛMac₁ + (1.0 - w) * θₛMac₂)
      
      Nsmooth = Nz - 1

      Kₛ_Smooth = zeros(Float64, Nsmooth)
      σ_Smooth = zeros(Float64, Nsmooth)
      Ψm_Smooth = zeros(Float64, Nsmooth)
      θs_Smooth = zeros(Float64, Nsmooth)
      θr_Smooth = zeros(Float64, Nsmooth)
      θsMac_Smooth = zeros(Float64, Nsmooth)

      for iZ=2:Nz
         Kₛ_Smooth[iZ-1] = Kₛ_AVER(hydro.Ks[iZ], hydro.Ks[iZ-1], 0.5)

         σ_Smooth[iZ-1] = σ_AVER(hydro.σ[iZ], hydro.σ[iZ-1], 0.5)

         Ψm_Smooth[iZ-1] = Ψₘ_AVER(hydro.Ψm[iZ], hydro.Ψm[iZ-1], 0.5)

         θs_Smooth[iZ-1] = θₛ_AVER(hydro.θs[iZ], hydro.θs[iZ-1], 0.5)

         θr_Smooth[iZ-1] = θᵣ_AVER(hydro.θr[iZ], hydro.θr[iZ-1], 0.5)

         θsMac_Smooth[iZ-1] = θₛMac_AVER(hydro.θsMac[iZ], hydro.θsMac[iZ-1], 0.5)
      end

   return
   end  # function: HYDRO_AVERAGE
   # ------------------------------------------------------------------
   
end  # module: hydroAver
# ............................................................