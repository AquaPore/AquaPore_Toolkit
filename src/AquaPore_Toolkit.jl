##========================================================================================
##                                                                                      ##
##                                 Soil Water ToolBox                                   ##
##                                                                                      ##
##========================================================================================

include("Including.jl")

module AquaPore_Toolkit
	import ..checking, ..hydrolabOpt, ..hydroStruct, ..infiltStart, ..ksModel, ..options, ..params, ..paths, ..plot, ..psdStart, ..reading, ..readSmap, ..rockFragment, ..smap2hypix, ..startKsModel, ..table, ..tableSmap, ..tool, ..wrc

	export AQUAPORE_TOOLBOX

	# ===============================================================
	#		FUNCTION : START_TOOLBOX
	# ==============================================================
	function AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="SoilWater", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="NewFormat")
		# _______________________ START: option/ param/ path _______________________ 

			Path_Home = @__DIR__
			Path_Data₀ = dirname(Path_Home)

			PathData_SoilWater = Path_Data₀ * "/data/INPUT/Data_SoilWater/" * SiteName_Soilwater

			if Soilwater_OR_Hypix⍰ == "SoilWater" # <>=<>=<>=<>=<>
				option = options.OPTIONS(PathData_SoilWater, SiteName_Soilwater)

				param = params.PARAM(PathData_SoilWater, SiteName_Soilwater)

				option.run.Hypix = false

			# elseif Soilwater_OR_Hypix⍰ == "Hypix"  # <>=<>=<>=<>=<>
			# 	hypixStart.HYPIX_START(SiteName_Soilwater)
			# 	printstyled("Run HyPix as indepent model")

			else  # <>=<>=<>=<>=<>
				error("Soilwater_OR_Hypix⍰ = $Soilwater_OR_Hypix⍰ not available needs to be either <SoilWater> or <Hypix>")
			end
				
			path = paths.PATH(option, PathData_SoilWater, SiteName_Soilwater)

		# ------------------------END: option/ param/ path---------------------------

		# ++++++++++++++++++++++ SCENARIOS ++++++++++++++++++++++++++++++++++++++++++
		N_Scenario = 1
		if option.run.Smap
			Scenarios = option.hydro.HydroModel_List
			N_Scenario =	length(Scenarios)
		end 
		for iSim =1:N_Scenario
			if option.run.Smap
				option.hydro.HydroModel⍰ = Scenarios[iSim]
				path = paths.PATH(option, PathData_SoilWater, SiteName_Soilwater)
				printstyled("\n +++++++++++++++++ SCENARIOS: option.hydro.HydroParam=$(option.hydro.HydroModel⍰)  $iSim / $N_Scenario, +++++++++++++++++ \n")
			end
		#..............................................................................


		# _______________________ START: reading _______________________ 
		printstyled("\n ----- START READING -----------------------------------------------\n")
	
			# DETERMINE WHICH SOILS/ PROFILE TO RUN: <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				IdSelect, IdSelect_True, Soilname, NiZ = reading.ID(PathIdSelect=path.inputSoilwater.IdSelect, PathOptionSelect=path.option.Select, PathModelName=path.option.ModelName)

				# Deriving opt parameters
					hydroₒ = hydroStruct.HYDROSTRUCT(option.hydro, 1)
					hydroₒ, optim = reading.HYDRO_PARAM(option.hydro, hydroₒ, 1, path.inputGuiSoilwater.GUI_HydroParam)

			# IF WE HAVE Θ(Ψ) DATA: <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				if option.data.θΨ && !(option.data.SimulationKosugiθΨK && option.hydro.HydroModel⍰ ≠"Kosugi" && option.hydro.σ_2_Ψm⍰=="Constrained")
					θ_θΨobs, Ψ_θΨobs, N_θΨobs = reading.θΨ(IdSelect, NiZ, path)

				elseif option.data.θΨ && option.data.SimulationKosugiθΨK && option.hydro.HydroModel⍰ ≠ "Kosugi" && option.hydro.σ_2_Ψm⍰=="Constrained" # Ading extra data
					try
						@info "\n	*** Reading θ(Ψ) data from $(path.tableSoilwater.TableComplete_θΨ) *** \n"
						θ_θΨobs, Ψ_θΨobs, N_θΨobs = reading.θΨ(IdSelect, NiZ, path)
					catch
						@warn "\n option.data.SimulationKosugiθΨK && option.hydro.HydroModel⍰ ≠:Kosugi && param.hydro.σ_2_Ψm⍰==Constrained => Kosugi simulation not performed yet! \n" 
						θ_θΨobs, Ψ_θΨobs, N_θΨobs = reading.θΨ(IdSelect, NiZ, path)
					end 		
				end  # if: option.data.θΨ


			# IF WE HAVE K(Θ) DATA: <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				if option.data.Kθ && !(option.data.SimulationKosugiθΨK && option.hydro.HydroModel⍰ ≠"Kosugi" && option.hydro.σ_2_Ψm⍰=="Constrained")
					K_KΨobs, Ψ_KΨobs, N_KΨobs = reading.KUNSATΨ(IdSelect, NiZ, path, path.inputSoilwater.Kunsat)

				elseif option.data.SimulationKosugiθΨK && option.hydro.HydroModel⍰ ≠ "Kosugi" 
					try
						@info "\n	*** Reading K(Ψ) data from $(path.tableSoilwater.Table_Smap_θΨK.csv) *** \n"
						println(path.tableSoilwater.TableComplete_KΨ)
						K_KΨobs, Ψ_KΨobs, N_KΨobs = reading.KUNSATΨ(IdSelect, NiZ, path, path.tableSoilwater.TableComplete_KΨ)
					catch
						@warn "\n *** option.data.SimulationKosugiθΨK && option.hydro.HydroModel⍰≠:Kosugi => Kosugi simulation not performed yet! *** \n"
						if "Ks" ∈ optim.ParamOpt
							K_KΨobs, Ψ_KΨobs, N_KΨobs = reading.KUNSATΨ(IdSelect, NiZ, path, path.inputSoilwater.Kunsat)
						end
					end # catch
				else
					K_KΨobs = []
					Ψ_KΨobs = []
					N_KΨobs = 1
				end  # if: Kθ			

			# IF WE HAVE THE HYDRAULIC PARAMETERS PRECOMPUTED FROM PREVIOUS SIMULATIONS	
				if option.run.HydroLabθΨ⍰ == "HydroParamPrecomputed"
					hydroₒ = hydroStruct.HYDROSTRUCT(option.hydro, 1)

					hydro, NiZ = tool.readWrite.READ_STRUCT_SIMPLE(hydroₒ, path.inputSoilwater.HydroParamPrecomputed)
					@info "\n	*** Reading hydro parameters from file *** \n "
				end # option.run.HydroLabθΨ⍰ == "HydroParamPrecomputed"


			# IF WE HAVE BULK DENSITY AND ROCK FRAGMENT DATA: <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				if option.data.Φ⍰ == "ρᵦ"
					RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil = reading.BULKDENSITY(IdSelect, NiZ, path.inputSoilwater.BulkDensity)

					# Φ  corrected for RockFragments
					Φ = rockFragment.ρᵦ_2_Φ(NiZ, option, RockFragment, ρₚ_Fine, ρₚ_Rock, ρᵦ_Soil)

				elseif option.data.Φ⍰ == "Φ" # Total Porosity
					RockFragment, Φ = reading.Φ(IdSelect, NiZ, path.inputSoilwater.Φ)
					
					Φ = rockFragment.injectRock.CORECTION_Φ!(NiZ, option, RockFragment, Φ)	
				end # option.data.Φ⍰ == :ρᵦ


			# IF WE HAVE INFILT DATA: <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				if option.data.Infiltration
					Tinfilt, ∑Infilt_Obs, N_Infilt, infiltParam = reading.INFILTRATION(IdSelect, NiZ, path.inputSoilwater.Infiltration, path.inputSoilwater.Infiltration_Param)
				end  # if: option.data.Infiltration


			# IF WE HAVE PSD DATA: <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				if option.data.Psd
					Rpart, ∑Psd, N_Psd = reading.PSD(IdSelect, NiZ, path.inputSoilwater.Psd)
				else
					∑Psd = []
				end  # if: option.data.Psd

				
			# IF WE WANT TO DERIVE Ks FROM θ(Ψ)
				if option.hydro.HydroModel⍰ == "Kosugi" && (option.run.KsModel || !(option.data.Kθ && "Ks" ∈ optim.ParamOpt))
					
					ksmodelτ, optimKsmodel = reading.KSΨMODEL_PARAM(NiZ, option, param, path.inputGuiSoilwater.GUI_KsModel) 
				end


			# IF WE HAVE PEDOLOGICAL⍰: <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				if option.data.Pedological⍰ == "Core"
					IsTopsoil, RockClass = reading.PEDOLOGICAL(IdSelect, NiZ, path.inputSoilwater.Pedological⍰)
				
				elseif option.data.Pedological⍰ == "Smap"
					IsTopsoil, Ks_Impermeable, RockClass, RockFragment, Smap_Depth, Smap_MaxRootingDepth, Smap_PermeabilityClass, Smap_SmapFH, Soilname = readSmap.SMAP(IdSelect_True, NiZ, path)
				end  # if: option.data.Pedological⍰


			#--- NON CORE ----
				# SMAP if we have information of the wetability of rocks:
					if option.data.RockWetability && option.run.Smap
						rfWetable = readSmap.ROCKFRAGMENT_WETTABLE(path.inputSmap.LookupTable_RockWetability)	
					end  # if: option.data.RockWetability

		printstyled("\n ----- END READING ----------------------------------------------- \n")
		
		# ------------------------END: reading---------------------------
		#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


		# _______________________ START: running HydroLabθΨ _______________________ 

		if option.run.HydroLabθΨ⍰ ≠ "No" && option.run.HydroLabθΨ⍰ ≠ "HydroParamPrecomputed"
		printstyled("\n ----- START RUNNING HYDROLABΘΨ ----------------------------------------------- \n"; color=:red)
			# STRUCTURES
				hydro = hydroStruct.HYDROSTRUCT(option.hydro, NiZ)
				hydroOther = hydroStruct.HYDRO_OTHERS(NiZ)
				hydro, optim = reading.HYDRO_PARAM(option.hydro, hydro, NiZ, path.inputGuiSoilwater.GUI_HydroParam; PrintScreen=true)

			# CHECKING THE DATA
				checking.CHECKING(option, option.hydro, optim)

			# TRANSFERING Φ -> hydro
				if option.data.Φ⍰ ≠ "No"
					for iZ =1:NiZ 
						hydro.Φ[iZ] = Φ[iZ]
					end
				end # option.data.Φ⍰ ≠ :No

			# CORRECT θ(Ψ) FOR ROCK FRAGMENT
			if option.run.RockCorection
				if option.rockFragment.RockInjectedIncluded⍰ =="InjectRock"
					@info "\n Correction for rock fragments  \n" 
					θ_θΨobs = rockFragment.injectRock.CORECTION_θΨ!(NiZ, N_θΨobs, RockFragment, θ_θΨobs)
				end #  option.rockFragment.RockInjectedIncluded⍰ ==:InjectRock

				if option.rockFragment.CorectStoneRockWetability
					@info "\n Correction for rock wettability  \n" 
					θ_θΨobs = rockFragment.CORECTION_θΨ_WETABLE!(NiZ, N_θΨobs, rfWetable, RockClass, RockFragment, θ_θΨobs, Ψ_θΨobs)
				end # option.rockFragment.CorrectStoneWetability
			end # if:option.run.RockCorection


			# OPTIMISING HYDRAULIC PARAMETERS
			if "Ks" ∈ optim.ParamOpt
				hydro, hydroOther = hydrolabOpt.HYDROLABOPT_START(NiZ=NiZ, ∑Psd=∑Psd, θ_θΨobs=θ_θΨobs, Ψ_θΨobs=Ψ_θΨobs, N_θΨobs=N_θΨobs, K_KΨobs=K_KΨobs, Ψ_KΨobs=Ψ_KΨobs, N_KΨobs=N_KΨobs, hydro=hydro, hydroOther=hydroOther, option=option, optionₘ=option.hydro, optim=optim, param=param)

			else
				hydro, hydroOther = hydrolabOpt.HYDROLABOPT_START(NiZ=NiZ, ∑Psd=∑Psd, θ_θΨobs=θ_θΨobs, Ψ_θΨobs=Ψ_θΨobs, N_θΨobs=N_θΨobs, hydro=hydro, hydroOther=hydroOther, option=option, optionₘ=option.hydro, optim=optim, param=param)

			end # "Ks" ∈ optim.ParamOpt

			# SPECIAL CASE
				if option.hydro.HydroModel⍰=="BrooksCorey" || option.hydro.HydroModel⍰=="ClappHornberger"
					for iZ=1:NiZ
						hydro.Ψga[iZ] = wrc.GREEN_AMPT(option.hydro, iZ, hydro)
					end
				end #  option.hydro.HydroModel⍰

		printstyled("\n ----- END: RUNNING HYDROLABΘΨ ----------------------------------------------- \n"; color=:green)
		end # option.run.HydroLabθΨ⍰
		# ------------------------END: running HydroLabθΨ--------------------------


		# _______________________ START: COMPUTE KS FROM Θ(Ψ) _______________________ 

			KₛModel = fill(NaN::Float64, NiZ)

			if option.hydro.HydroModel⍰ == "Kosugi" && (option.run.KsModel || !(option.data.Kθ && "Ks" ∈ optim.ParamOpt))
				printstyled("\n ----- START RUNNING Ks Model from θ(Ψ)  ----------------------------------------------- \n"; color=:red)
					printstyled("		Running KsModel= ", option.ksModel.KₛModel⍰, "\n" ; color=:green)

				# Default value
					🎏_IsTopsoil=false; 🎏_RockFragment=false; IsTopsoil₀=[]; RockFragment₀=[]; Ks_Impermeable₀=[]

				if  @isdefined RockFragment
               🎏_RockFragment = true
               RockFragment₀   = RockFragment
				end
				
				if @isdefined IsTopsoil
               🎏_IsTopsoil = true
               IsTopsoil₀   = IsTopsoil
				end

				if option.run.Smap
               🎏_IsTopsoil    = true
               🎏_RockFragment = true
               IsTopsoil₀      = IsTopsoil
               RockFragment₀   = RockFragment
               Ks_Impermeable₀ = Ks_Impermeable
				end

				hydro, KₛModel, N_Class = startKsModel.START_KSΨMODEL(hydro, KₛModel, ksmodelτ, NiZ, optim, optimKsmodel, option, param, path; 🎏_IsTopsoil=🎏_IsTopsoil, 🎏_RockFragment=🎏_RockFragment, IsTopsoil=IsTopsoil₀, RockFragment=RockFragment₀, Ks_Impermeable=Ks_Impermeable₀, ∑Psd=∑Psd)

				printstyled("\n ----- END RUNNING Ks Modelfrom θ(Ψ) ----------------------------------------------- \n";color=:green)
			end # if: option.hydro.HydroModel⍰ == :Kosugi
		# ------------------------END:  COMPUTE KS FROM Θ(Ψ) -------------------------- 


		# _______________________ START: IntergranularMixingPsd _______________________ 
		if option.run.IntergranularMixingPsd 
			printstyled("\n ----- START RUNNING IntergranularMixingPsd ----------------------------------------------- \n";color=:green)
			# STRUCTURES
				hydroPsd = hydroStruct.HYDROSTRUCT(option.psd, NiZ)
				hydroOther_Psd = hydroStruct.HYDRO_OTHERS(NiZ)
				hydroPsd, optim_Psd = reading.HYDRO_PARAM(option.psd, hydroPsd, NiZ, path.inputGuiSoilwater.GUI_HydroParam)

			# CHECKING THE DATA
				checking.CHECKING(option, option.psd, optim)

			# TRANSFERING Φ -> hydro
				for iZ =1:NiZ 
					hydroPsd.Φ[iZ] = Φ[iZ]
				end

			# PSD model
			if @isdefined hydro
				paramPsd, N_Psd, θ_Rpart, Ψ_Rpart, Psd, hydroPsd = psdStart.START_PSD(∑Psd=∑Psd, hydro=hydro, hydroPsd=hydroPsd, NiZ=NiZ, N_Psd=N_Psd, N_θΨobs=N_θΨobs, option=option, param=param, Rpart=Rpart, θ_θΨobs=θ_θΨobs, Ψ_θΨobs=Ψ_θΨobs)
			else
				paramPsd, N_Psd, θ_Rpart, Ψ_Rpart, Psd, hydroPsd = psdStart.START_PSD(∑Psd=∑Psd, hydroPsd=hydroPsd, NiZ=NiZ, N_Psd=N_Psd, N_θΨobs=N_θΨobs, option=option, param=param, Rpart=Rpart, θ_θΨobs=θ_θΨobs, Ψ_θΨobs=Ψ_θΨobs)
			end

			# Deriving the hydraulic parameters of PSD
				KunsatModel_Psd = fill(0.0::Float64, NiZ)

				hydroPsd, hydroOther_Psd = hydrolabOpt.HYDROLABOPT_START(NiZ=NiZ, ∑Psd=∑Psd, θ_θΨobs=θ_Rpart, Ψ_θΨobs=Ψ_Rpart, N_θΨobs=N_Psd, hydro=hydroPsd, hydroOther=hydroOther_Psd, option=option, optionₘ=option.psd, optim=optim_Psd, param=param) 

				printstyled("\n 	----- START RUNNING Ks Model from θ(Ψ)PSD  -----------------------------------------------"; color=:green)
					if  (@isdefined RockFragment) && (@isdefined IsTopsoil)
						hydroPsd, KₛModel = startKsModel.START_KSΨMODEL(hydroPsd, option, param, path, KₛModel, path.option.ModelName, ksmodelτ, NiZ, optim, optimKsmodel; 🎏_IsTopsoil=true, 🎏_RockFragment=true, IsTopsoil=IsTopsoil, RockFragment=RockFragment)
					
					elseif (@isdefined RockFragment) && !(@isdefined IsTopsoil)	
						hydroPsd, KₛModel = startKsModel.START_KSΨMODEL(hydroPsd, option, param, path, KₛModel, path.option.ModelName, ksmodelτ, NiZ, optim, optimKsmodel; 🎏_RockFragment=true, RockFragment=RockFragment)
					
					elseif !(@isdefined RockFragment) && (@isdefined IsTopsoil)
						hydroPsd, KₛModel = startKsModel.START_KSΨMODEL(hydroPsd, option, param, path, KₛModel, path.option.ModelName, ksmodelτ, NiZ, optim, optimKsmodel; 🎏_IsTopsoil=true, IsTopsoil=IsTopsoil)
					
					elseif !(@isdefined RockFragment) && !(@isdefined IsTopsoil)
						hydroPsd, KₛModel = startKsModel.START_KSΨMODEL(hydroPsd, option, param, path, KₛModel, path.option.ModelName, ksmodelτ, NiZ, optim, optimKsmodel)
					end # if: RockFragment && IsTopsoil
				printstyled("\n 	----- END RUNNING Ks Model from θ(Ψ)PSD  ----------------------------------------------- \n"; color=:green)					
			
			printstyled("\n ----- END: RUNNING IntergranularMixingPsd ----------------------------------------------- \n"; color=:yellow)
		end
		# ------------------------END: IntergranularMixingPsd---------------------------  


		# _______________________ START: Infiltration _______________________ 

		if option.run.Infiltration
			printstyled("\n ----- START RUNNING INFILTRATION -----------------------------------------------"; color=:green)
			# STRUCTURES
				hydroInfilt = hydroStruct.HYDROSTRUCT(option.infilt, NiZ)
				hydroOther_Infilt = hydroStruct.HYDRO_OTHERS(NiZ)
				hydroInfilt, optim_Infilt = reading.HYDRO_PARAM(option.psd, hydroInfilt, NiZ, path.inputGuiSoilwater.GUI_HydroParam)

			# CHECKING THE DATA
				checking.CHECKING(option, option.infilt, optim)

			# TRANSFERING Φ -> hydro
				for iZ =1:NiZ 
					hydroInfilt.Φ[iZ] = Φ[iZ]
				end

			# RUNNING INFILTRATION MODEL
			if @isdefined hydro
				infiltOutput, hydroInfilt, ∑Infilt_3D, ∑Infilt_1D = infiltStart.START_INFILTRATION(∑Infilt_Obs=∑Infilt_Obs, ∑Psd=∑Psd, hydro=hydro, hydroInfilt=hydroInfilt, infiltParam=infiltParam, N_Infilt=N_Infilt, NiZ=NiZ, option=option, param=param,Tinfilt=Tinfilt)
			else
				infiltOutput, hydroInfilt, ∑Infilt_3D, ∑Infilt_1D = infiltStart.START_INFILTRATION(∑Infilt_Obs=∑Infilt_Obs, ∑Psd=∑Psd, hydroInfilt=hydroInfilt, infiltParam=infiltParam, N_Infilt=N_Infilt, NiZ=NiZ, option=option, param=param, Tinfilt=Tinfilt)
			end
		printstyled("\n ----- END: RUNNING Infiltration ----------------------------------------------- \n"; color=:yellow)
		end # option.run.Infiltration

		# ------------------------END: Infiltration---------------------------

			
		# _______________________ START: Smap_2_HyPix ______________________
		if option.run.Smap2Hypix
		printstyled("\n----- START RUNNING Smap ➡ Smap_Hypix  ----------------------------------------------\n"; color=:green)
			smap2hypix.SMAP_2_HYPIX(hydro, NiZ, option.hydro, param, path, Smap_Depth, Smap_MaxRootingDepth, Soilname)
		printstyled("\n----- END RUNNING  Smap ➡ Smap_Hypix  ----------------------------------------------- , \n"; color=:yellow)
		end  # if: Smap2Hypix 

		# ------------------------END: Smap_2_HyPix---------------------------
		
		# _______________________ START: Temporary _______________________ 

		# if option.run.Temporary
		# 	temporary.KS_SMAP()
		# end

		# ------------------------END: Temporary---------------------------  


		#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# _______________________ START: table _______________________ 

		# if option.run.ChangeHydroModel
		# 	table.hydroLab.TABLE_EXTRAPOINTS_Kθ(option.hydro, hydro, IdSelect, param.hydro.K_Table, NiZ, path.inputSoilwater.Kunsat)
			
		# 	table.hydroLab.TABLE_EXTRAPOINTS_θΨ(option.hydro, hydro, IdSelect, NiZ, path.inputSoilwater.Ψθ , param.hydro.TableComplete_θΨ; Orientation="Vertical")
		# end

		if option.run.HydroLabθΨ⍰ ≠ "No" && option.run.HydroLabθΨ⍰ ≠ "HydroParamPrecomputed" # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
			# CORE OUTPUT
				table.hydroLab.θΨK(hydro, hydroOther, IdSelect[1:NiZ], KₛModel[1:NiZ], NiZ, path.tableSoilwater.Table_θΨK)

				# When optimising other model than Kosugi we do not have a model for σ_2_Ψm⍰. Therefore we assume that θ(Ψ) and K(θ) derived by Kosugi from very dry to very wet are physical points
				# if option.hydro.HydroModel⍰ == "Kosugi" && option.hydro.σ_2_Ψm⍰=="Constrained"
				if option.run.ChangeHydroModel
				println(path.tableSoilwater.TableComplete_KΨ)
					table.hydroLab.TABLE_EXTRAPOINTS_Kθ(option.hydro, hydro, IdSelect, param.hydro.K_Table, NiZ, path.tableSoilwater.TableComplete_KΨ)
			
					table.hydroLab.TABLE_EXTRAPOINTS_θΨ(option.hydro, hydro, IdSelect, NiZ, path.tableSoilwater.TableComplete_θΨ, param.hydro.TableComplete_θΨ; Orientation="Vertical")
				end # option.run.ChangeHydroModel
				# end # if: option.hydro.HydroModel⍰ == :Kosugi && option.hydro.σ_2_Ψm⍰ == Constrained

				# IF SMAP OUTPUTS
				if option.run.Smap
					tableSmap.θΨK(hydro, hydroOther, IdSelect, KₛModel, NiZ, path.tableSmap.Table_θΨK, Smap_Depth, Soilname)

					# When all the models are performed
					if iSim==length(Scenarios)
						tableSmap.SMAP(hydro, IdSelect, IsTopsoil, NiZ, option.hydro, param, path, RockFragment, Smap_Depth, Smap_MaxRootingDepth, Smap_PermeabilityClass, Smap_SmapFH, Soilname)
					end
				end # option.run.Smap	
			end # option.run.HydroLabθΨ⍰ ≠ :No && option.run.HydroLabθΨ⍰ ≠ :File

			if option.run.KsModel # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				table.ksmodel.KSMODEL(hydro, IdSelect, KₛModel,  path.tableSoilwater.Table_KsModel)
				table.ksmodel.KSMODEL_τ(ksmodelτ, N_Class, path.tableSoilwater.Table_KsModel_τ)
			end  # if: option.run.KsModel

			if option.run.Infiltration # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				table.infilt.HYDRO_INFILT(hydroInfilt, IdSelect, NiZ, path.tableSoilwater.Table_HydroInfilt)

				table.infilt.INFILT(IdSelect, NiZ, infiltOutput, path.tableSoilwater.Table_Infilt)
			end # option.run.Infiltration


			if option.run.IntergranularMixingPsd # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
				table.psd.PSD(IdSelect[1:NiZ], NiZ, paramPsd, path.tableSoilwater.Table_Psd)

				table.psd.θΨK_PSD(hydroPsd, IdSelect, KunsatModel_Psd, NiZ, path.tableSoilwater.Table_Psd)
				
				if option.psd.Table_Psd_θΨ_θ
					table.psd.PSD_θΨ_θ(IdSelect, hydroPsd, NiZ, option, param, path.tableSoilwater.Table_Psd_θΨ_θ)
				end
			end # option.run.IntergranularMixingPsd

		# ------------------------END: table---------------------------
		#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
		
		# _______________________ START: plotting _______________________ 

			printstyled("\n		=== START: PLOTTING  === \n";color=:green)

			# Checking the maximum number of plotting
				param.globalparam.N_iZ_Plot_End = min(param.globalparam.N_iZ_Plot_End, NiZ)

				if option.run.HydroLabθΨ⍰ ≠ "No" && option.run.HydroLabθΨ⍰ ≠ "HydroParamPrecomputed" && option.hydro.Plot_θΨ
					plot.lab.HYDROPARAM(hydro, hydroOther, IdSelect, K_KΨobs, NiZ, N_KΨobs, N_θΨobs, optim, option, param, path, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs)
				end

				if option.run.IntergranularMixingPsd
					if option.psd.Plot_θr && option.run.HydroLabθΨ⍰ ≠ "No" # <>=<>=<>=<>=<>
						plot.psd.PLOT_θr(∑Psd, hydro, hydroPsd, NiZ, param, path.plotSoilwater.Plot_Psd_θr)
					end
					if option.psd.Plot_IMP_Model # <>=<>=<>=<>=<>
						plot.psd.PLOT_IMP_MODEL(∑Psd, hydro, IdSelect, NiZ, N_Psd, option, param, path.plotSoilwater.Plot_IMP_model, Psd, Rpart) 
					end
					if option.psd.Plot_Psd_θΨ # <>=<>=<>=<>=<>
						plot.psd.PLOT_PSD_θΨ(hydro, hydroPsd, IdSelect, NiZ, N_Psd, N_θΨobs, option, param, path.plotSoilwater.Plot_Psd_θΨ, θ_Rpart, θ_θΨobs, Ψ_Rpart, Ψ_θΨobs)
					end
				end # option.run.IntergranularMixingPsd

				if option.run.Infiltration # <>=<>=<>=<>=<>
					if option.infilt.Plot_∑Infiltration  
						plot.infilt.PLOT_∑INFILT(∑Infilt_1D, ∑Infilt_3D, ∑Infilt_Obs, IdSelect, infiltOutput, N_Infilt, NiZ, option, param, path.plotSoilwater.Plot_∑infilt_Opt, Tinfilt)
					end
					if option.infilt.Plot_θΨ
						if option.run.HydroLabθΨ⍰ ≠ "No"
							plot.infilt.PLOT_∑INFILT_θΨ(hydroInfilt, IdSelect, NiZ, optim, option, param, path.plotSoilwater.Plot_∑infilt_θΨ; hydro=hydro)
						else
							plot.infilt.PLOT_∑INFILT_θΨ(hydroInfilt, IdSelect, NiZ, optim, option, param, path.plotSoilwater.Plot_∑infilt_θΨ)
						end # option.run.HydroLabθΨ⍰
					end # option.run.Infiltration
				end # option.run.Infiltration

				# if option.run.Smap # <>=<>=<>=<>=<>
				# 	plotSmap.makie.HYDROPARAM(hydro, IdSelect, K_KΨobs, KₛModel, NiZ, N_KΨobs, N_θΨobs, option, path, Smap_Depth, Soilname, θ_θΨobs, Ψ_KΨobs, Ψ_θΨobs)
				# end # option.run.Smap
			
			printstyled("\n		=== END: PLOTTING  === \n"; color=:green)
		
		# ------------------------END: plotting---------------------------  
		#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

		end #iSim

	end  # function: START_TOOLBOX
	# ..............................................................

end # module soilwater_toolbox # module soilwater_toolbox

printstyled("\n\n ===== START SOIL WATER TOOLBOX =====, \n"; color=:green)
	
	# @time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="SoilWater", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="NewFormat")
	
	# @time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="SoilWater", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="Nsdr")

	# @time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="SoilWater", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="Int")


	# @time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="SoilWater", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="SFF")

	# @time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="SoilWater", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="SmapSwat")

	# @time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="SoilWater", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="TestSmapHydro20220728")

	# @time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="SoilWater", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="Lysimeters")

	#  @time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="Hypix", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="Convert")

	# @time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="Hypix", SiteName_Hypix="TESTCASE", SiteName_Soilwater="Convert")
	# @time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="SoilWater", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="SmapSmapNZSnapshot20210823")
	
	# @time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="SoilWater", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="Unsoda")

	@time AquaPore_Toolkit.AQUAPORE_TOOLBOX(;Soilwater_OR_Hypix⍰="SoilWater", SiteName_Hypix="LYSIMETERS", SiteName_Soilwater="SmapHydro")

printstyled("\n ==== END SOIL WATER TOOLBOX ====, \n"; color=:red)