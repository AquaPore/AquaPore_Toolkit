# =============================================================
#		module: θaver
# =============================================================
module θaver
   import ..discretisation
   export θAVER

   # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #		FUNCTION : θ_AVERAGE
   # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #    function θAVER(discret; Z=Z, θ_Reduced=θ_Reduced, Nz=Nz, Nit_Reduced=Nit_Reduced, Zaver=400.0)
   #       # Make sure that Zaver is whitin physical bounds
   #          Zaver = min(Zaver, Z[Nz])

   #       # Memory
   #          θsim_Aver = fill(0.0::Float64, Nit_Reduced)

   #       # For every time step
   #       for iT = 1:Nit_Reduced
   #          iZ_Max = 1
   #          for iZ = 1:Nz
   #             if Z[iZ] ≤ Zaver
   #                θsim_Aver[iT] += θ_Reduced[iT, iZ] * discret.ΔZ[iZ]
   #             else
   #                iZ_Max = iZ
   #                break
   #             end
   #             iZ_Max = iZ
   #          end # iZ=1:Nz
   #          if Z[iZ_Max-1] + eps(100.0) < Zaver
   #             θsim_Aver[iT] += θ_Reduced[iT, iZ_Max] * (Zaver - Z[iZ_Max-1])
   #          end
            
   #          θsim_Aver[iT] =  θsim_Aver[iT] / Zaver
   #       end # for iT

   #    return θsim_Aver
   #    end  # function: θ_AVERAGE
   # # ...................................................................

   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θ_AVERAGE2
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function θAVER(discret; Z=Z, θ_Reduced=θ_Reduced, Nz=Nz, Nit_Reduced=Nit_Reduced, ZaverArray=[300.0, 600.0, 1000.0])
         # Make sure that Zaver is whitin physical bounds


         # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            function STORAGE_0_Z(discret; Z=Z, θ_Reduced=θ_Reduced, Nz=Nz, Nit_Reduced=Nit_Reduced, Zaver=Zaver)
               # Make sure that Zaver is whitin physical bounds
               Zaver = min(Zaver, Z[Nz])

               # Memory
                  Storage_1D = fill(0.0::Float64, Nit_Reduced)
      
               # For every time step
               for iT = 1:Nit_Reduced
                  iZ_Max = 1
                  for iZ = 1:Nz
                     if Z[iZ] ≤ Zaver
                        Storage_1D[iT] += θ_Reduced[iT, iZ] * discret.ΔZ[iZ]
                     else
                        iZ_Max = iZ
                        break
                     end
                     iZ_Max = iZ
                  end # iZ=1:Nz
                  if Z[iZ_Max-1] + eps(100.0) < Zaver
                     Storage_1D[iT] += θ_Reduced[iT, iZ_Max] * (Zaver - Z[iZ_Max-1])
                  end
                  
                  # Storage_1D[iT] =  Storage_1D[iT] / Zaver
               end # for iT
                  
            return Storage_1D
            end  # function: \theta averag 
         # ------------------------------------------------------------------

         Nzaverage_Array = length(ZaverArray)

         Storage_Aver = fill(0.0::Float64, Nit_Reduced, Nzaverage_Array)
         θsim_Aver    = fill(0.0::Float64, Nit_Reduced, Nzaverage_Array)

         for iZaver =1:Nzaverage_Array
            Storage_Aver[:, iZaver] = STORAGE_0_Z(discret; Z=Z, θ_Reduced=θ_Reduced, Nz=Nz, Nit_Reduced=Nit_Reduced, Zaver=ZaverArray[iZaver])
         end

         θsim_Aver[:, 1] =  Storage_Aver[:, 1] / ZaverArray[1]
         if Nzaverage_Array ≥ 2
            for iZaver = Nzaverage_Array:-1:2
               θsim_Aver[:,iZaver] = (Storage_Aver[:,iZaver] - Storage_Aver[:,iZaver-1]) / (ZaverArray[iZaver] - ZaverArray[iZaver-1])
            end
         end
      return θsim_Aver
      end  # function: θ_AVERAGE
   # ----------------------------------------------------------------

   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : ∑QZ
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function ∑QZ(∑Q_Z, iOpt_Count, Nz, Z, ΔQ_Reduced; Zq=600.0) 
         iZ = 1
         while iZ ≤ Nz && Z[iZ] ≤ Zq 
            iZ +=1
         end
         N_iRoot = min(iZ-1, Nz)
         
         ∑Q_Z[iOpt_Count] = sum(ΔQ_Reduced[:, N_iRoot])
      return ∑Q_Z
      end  # function: Q600
      # ------------------------------------------------------------------

end  # module: θaver
# ............................................................

