 # =============================================================
#		MODULE: optimize
# =============================================================
 module optimizeOptim
   export SEARCHRANGE_OPTIM
   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SEARCHRANGE
	#		Required by OPTIM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function SEARCHRANGE_OPTIM(optionₘ, optim)
         ParamOpt_Min₂ = copy(optim.ParamOpt_Min)
         ParamOpt_Max₂ = copy(optim.ParamOpt_Max)

         # Log transform
         for iParam=1:optim.NparamOpt
            if optim.ParamOpt_LogTransform[iParam]
               ParamOpt_Min₂[iParam] = log1p(optim.ParamOpt_Min[iParam])
               ParamOpt_Max₂[iParam] = log1p(optim.ParamOpt_Max[iParam])
            end
         end

      # Making sure that for constrained optimisation Ψm is between 0 & 1
         if (optionₘ.σ_2_Ψm⍰=="Constrained") && ("Ψm" ∈ optim.ParamOpt)
            iψm = findfirst(isequal("Ψm"), optim.ParamOpt)[1]

            ParamOpt_Min₂[iψm] = 0.0
            ParamOpt_Max₂[iψm] = 1.0
         end

         Lower = Float64.(ParamOpt_Min₂)
         Upper = Float64.(ParamOpt_Max₂)

      return Lower, Upper
      end  # function: SEARCHRANGE
end # module optimize