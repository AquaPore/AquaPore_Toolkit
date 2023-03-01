# =============================================================
#		module option
# ===========================================================
module options

	using Configurations, TOML

	# What available data we have?
	@option struct DATA
		Infiltration::Bool
		Kθ::Bool
		Pedological⍰::String
		Psd::Bool
		RockWetability::Bool
		SimulationKosugiθΨK::Bool
		θΨ::Bool
		Φ⍰::String
	end # struct DATA

	# What model wanting to run
	@option mutable struct RUN
      ChangeHydroModel::Bool
      HydroLabθΨ⍰::String
      Hypix::Bool
      Infiltration::Bool
      IntergranularMixingPsd ::Bool
		KsModel::Bool
      RockCorection::Bool
		Smap::Bool
		Smap2Hypix::Bool
      Temporary::Bool
	end

	@option mutable struct OTHER
		PlotVscode::Bool
	end
	@option mutable struct SMAP
		Nothings
	end
	@option mutable struct HYDRO
		HydroModel⍰::String
		HydroModel_List
		θrOpt⍰::String
		σ_2_Ψm⍰::String
		Plot_θΨ::Bool
	end

	@option mutable struct KSMODEL
      Class                 ::Bool
      KₛModel⍰              ::String
      OptIndivSoil          ::Bool
      Opt_Kθ                ::Bool
      Plot_KsModel          ::Bool
	end

	@option mutable struct PSD
		Model⍰::String
		OptimizePsd⍰::String
		Psd_2_θr⍰::String
		∑Psd_2_ξ1::Bool
		HydroModel⍰::String
		θrOpt⍰::String
		σ_2_Ψm⍰::String
		Plot_Psd_θΨ::Bool
		Plot_θr::Bool
		Plot_IMP_Model::Bool
		Table_Psd_θΨ_θ::Bool
	end
	@option mutable struct ROCKFRAGMENT
		CorectStoneRockWetability::Bool
		RockInjectedIncluded⍰::String
	end
	@option mutable struct INFILT
		DataSingleDoubleRing⍰::String  
		OptimizeRun⍰::String  
		Model⍰::String
		BestUniv_Continous::Bool                	  
		SorptivitySplitModel⍰::String  
		SorptivityModel⍰::String     
		HydroModel⍰::String             
		θrOpt⍰::String            
		σ_2_Ψm⍰::String                   	 
		Plot_∑Infiltration::Bool           
		Plot_θΨ::Bool                
	end

		@option mutable struct OPTION
			data::DATA
			hydro::HYDRO
			infilt::INFILT
			general::OTHER
			psd::PSD
			rockFragment::ROCKFRAGMENT
			run::RUN
			smap::SMAP
			ksModel::KSMODEL
		end # struct OPTION
	
	#__________________________________________________________________
	#..................................................................

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPTION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OPTIONS(Path_Data, SiteName)	
			Path = Path_Data * "/ParamOptionPath/" * SiteName * "_Option.toml"
			@assert isfile(Path)
		return Configurations.from_toml(OPTION, Path)
		end  # function: OPTION

end # module option 
# end OPTION