# =============================================================
#		module: checking
# =============================================================
module checking

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : CHECKING
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function CHECKING(option, optionₘ, optim)


		 # ------------ Cannot run simultaneously HyIx and soilwater---------------------
		 	if option.run.Hypix && (option.run.ChangeHydroModel || option.run.HydroLabθΨ⍰≠"No" || option.run.Infiltration || option.run.IntergranularMixingPsd ||  option.run.Smap2Hypix || option.run.Temporary)
				error("*** Cannot run simulataneously HyPix and SoilWater ***")

        # ------------ CHECKING HydroLabθΨ---------------------
			elseif "Ks" ∉ optim.ParamOpt && option.data.Kθ && option.run.HydroLabθΨ⍰ == "Opt" 
				@warn("*** If  optim.ParamOpt && option.data.Kθ && option.run.HydroLabθΨ⍰==Opt ⇒ Ks ∈ optim.ParamOpt  ***")

		 	elseif "Ks" ∈ optim.ParamOpt && !(option.data.Kθ) 
				error("*** If Ks ∈ optim.ParamOpt ⇒ option.data.Kθ ***")

			elseif option.run.HydroLabθΨ⍰≠"No" && !option.data.θΨ
				error("*** If option.run.HydroLabθΨ⍰ ⇒ option.data.θΨ ***")
			
			elseif option.hydro.ΨmacMat_2_σMac_ΨmMac && ("ΨmMac" ∈ optim.ParamOpt || "σMac" ∈ optim.ParamOpt)
				error("*** If option..hydro.ΨmacMat_2_σMac_ΨmMac ⇒ ΨmMac || σMac param should not be optimised pls change in GUI_HydroParam.csv ***")

			elseif optionₘ.HydroModel⍰=="Kosugi" && "θsMacMat" ∈ optim.ParamOpt
				error("*** optionₘ.HydroModel⍰==Kosugi && optionₘ.HydroModel⍰==Bimodal THAN optionₘ.HydroModel⍰ == Φ ***")
							              
			elseif optionₘ.σ_2_Ψm⍰ =="UniqueRelationship" && "Ψm" ∈ optim.ParamOpt
				error("*** optionₘ.σ_2_Ψm⍰ ==UniqueRelationship THAN Ψm does not need to be optimised ***")
			
			elseif optionₘ.HydroModel⍰ == "Kosugi" && (optionₘ.σ_2_Ψm⍰ =="Constrained" && "Ψm" ∉ optim.ParamOpt) 
				error("*** optionₘ.σ_2_Ψm⍰ ==Constrained THAN Ψm needs to be optimised ***")

			elseif  (optionₘ.θrOpt⍰=="σ_2_θr") && ("θr" ∈ optim.ParamOpt)
				error("*** optionₘ.θrOpt⍰==σ_2_θr THAN θr does not need to be optimized ***") 

			elseif (optionₘ.θrOpt⍰=="σ_2_θr") && ("σ" ∉ optim.ParamOpt)
				error("*** optionₘ.θrOpt⍰==σ_2_θr THAN σ needs to be optimized ***")

			elseif  (optionₘ.θrOpt⍰=="ParamPsd") && ("θr"∉ optim.ParamOpt) && !(option.data.Psd) # Derive θr frpm PSD
				error("*** optionₘ.θrOpt⍰==ParamPsd THAN option.run.IntergranularMixingPsd=true ***")

        	elseif option.data.SimulationKosugiθΨK && "Ks" ∉ optim.ParamOpt
            error("***  Ks  ∉ optim.ParamOpt && option.smap.SimulationKosugiθΨK THAN Ks ∈ optim.ParamOpt ***")

			elseif option.run.HydroLabθΨ⍰ == "HydroParamPrecomputed" 
				@warn("***  Running from HydroParamPrecomputed.csv ***")

			elseif option.hydro.ΨmacMat_2_σMac_ΨmMac && ("ΨmMac" ∈ optim.ParamOpt || "σMac" ∈ optim.ParamOpt )
				error("***  ΨmacMat_2_σMac_ΨmMac = true, ΨmMac ∈ optim.ParamOpt || σMac ∈ optim.ParamOpt ***") 

		# ------------ CHECKING SimulationKosugiθΨK which simulates other θ(ψ) & K(ψ) functions from Kosugi  θ(ψ) & K(ψ) models --------------------
			elseif option.data.SimulationKosugiθΨK && option.run.RockCorection
				error("***  option.data.SimulationKosugiθΨK && option.run.RockCorection THAN rock fragments are double corrected ***")

		# ------------ CHECKING Infiltration model--------------------
			elseif option.run.Infiltration && !(option.data.Infiltration)
				error("***  option.run.Infiltration ⇒option.data.Infiltration ***")

		# ------------ CHECKING Particle Size Distribution model--------------------
			elseif option.run.IntergranularMixingPsd && !(option.data.Psd)
				error("***  option.run.IntergranularMixingPsd ⇒option.data.Psd ***")

			elseif option.run.IntergranularMixingPsd && option.data.Φ⍰=="No"
				error("***  option.run.IntergranularMixingPsd ⇒ option.data.Φ⍰ ≠ No ***")

			elseif option.run.IntergranularMixingPsd && "Ks" ∈ optim.ParamOpt
					error("*** option.run.IntergranularMixingPsd ⇒ Ks ∉ optim.ParamOpt ***")
			
				# ------------ CHECKING Smap--------------------
			elseif option.run.Smap && option.data.Pedological⍰≠"Smap"
				error("*** If option.run.Smap ⇒option.data.Pedological==Smap ***")

			elseif option.rockFragment.CorectStoneRockWetability && !(option.run.Smap)
				@warn("*** If option.data.RockWetability ⇒ option.run.Smap ***")
			
		# ------------ CHECKING KsModel---------------------
			elseif option.run.KsModel && !(option.data.θΨ)
				error("*** If option.run.KsModel ⇒ option.data.θΨ == true ***")

		# ------------ CHECKING Smap_2_Hypix---------------------
			elseif option.run.Smap2Hypix && (option.data.Pedological⍰≠"Smap")
				error("*** option.run.Smap2Hypix ⇒option.data.Pedological⍰ == Smap ***")

			elseif option.run.Smap2Hypix && !(option.run.Smap) 
				error("*** option.run.Smap2Hypix ⇒ option.run.Smap  ➡ true***")

			elseif option.rockFragment.CorectStoneRockWetability && !(option.data.RockWetability)
				error("*** option.data.RockWetability ⇒ option.rockFragment.CorectStoneRockWetability  ➡ true***")

		end # Check error
    	return nothing
    	end  # function: CHECKING
   
end  # module: checking
# ............................................................