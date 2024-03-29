# =============================================================
#		MODULE: hydroStruct
# =============================================================
module hydroStruct
	import ..tool
	export HYDROSTRUCT
   using Base

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : HYDRAULIC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Base.@kwdef mutable struct KOSUGI  # <>=<>=<>=<>=<>=<>=<>=<>=<>
         θs             :: Vector{Float64}
         θr             :: Vector{Float64}
         Ks             :: Vector{Float64}
         Ψm             :: Vector{Float64}
         σ              :: Vector{Float64}
         θsMacMat_ƞ     :: Vector{Float64}
         σMac           :: Vector{Float64}
         ΨmMac          :: Vector{Float64}
         ΨmacMat        :: Vector{Float64}
         So             :: Vector{Float64}
         θsMacMat       :: Vector{Float64}
         Φ              :: Vector{Float64}
         τa             :: Vector{Float64}
         τb             :: Vector{Float64}
         τc             :: Vector{Float64}
         τₚ             :: Vector{Float64}
         τaMac          :: Vector{Float64}
         τbMac          :: Vector{Float64}
         τcMac          :: Vector{Float64}
			
         θs_Min         :: Vector{Float64}
         θr_Min         :: Vector{Float64}
         Ks_Min         :: Vector{Float64}
         Ψm_Min         :: Vector{Float64}
         σ_Min          :: Vector{Float64}
         θsMacMat_ƞ_Min :: Vector{Float64}
         σMac_Min       :: Vector{Float64}
         ΨmMac_Min      :: Vector{Float64}
         ΨmacMat_Min    :: Vector{Float64}
         So_Min         :: Vector{Float64}
         θsMacMat_Min   :: Vector{Float64}
         Φ_Min          :: Vector{Float64}
         τa_Min         :: Vector{Float64}
         τb_Min         :: Vector{Float64}
         τc_Min         :: Vector{Float64}
         τₚ_Min         :: Vector{Float64}
         τaMac_Min      :: Vector{Float64}
         τbMac_Min      :: Vector{Float64}
         τcMac_Min      :: Vector{Float64}

         θs_Max         :: Vector{Float64}
         θr_Max         :: Vector{Float64}
         Ks_Max         :: Vector{Float64}
         Ψm_Max         :: Vector{Float64}
         σ_Max          :: Vector{Float64}
         θsMacMat_ƞ_Max :: Vector{Float64}
         σMac_Max       :: Vector{Float64}
         ΨmMac_Max      :: Vector{Float64}
         ΨmacMat_Max    :: Vector{Float64}
         So_Max         :: Vector{Float64}
         θsMacMat_Max   :: Vector{Float64}
         Φ_Max          :: Vector{Float64}
         τa_Max         :: Vector{Float64}
         τb_Max         :: Vector{Float64}
         τc_Max         :: Vector{Float64}
         τₚ_Max         :: Vector{Float64}
         τaMac_Max      :: Vector{Float64}
         τbMac_Max      :: Vector{Float64}
         τcMac_Max      :: Vector{Float64}
		end # struct KOSUGI

		
		Base.@kwdef mutable struct VANGENUCHTEN # <>=<>=<>=<>=<>=<>=<>=<>=<>
         θs          ::	Vector{Float64}
         θr          ::	Vector{Float64}
         N           ::	Vector{Float64}
         Ψvg         ::	Vector{Float64}
         ΨmacMat     :: Vector{Float64}
         Ks          ::	Vector{Float64}
         Km          ::	Vector{Float64}
         Φ           ::	Vector{Float64}
			
         θs_Min      ::	Vector{Float64}
         θr_Min      ::	Vector{Float64}
         N_Min       ::	Vector{Float64}
         Ψvg_Min     ::	Vector{Float64}
         ΨmacMat_Min :: Vector{Float64}
         Ks_Min      ::	Vector{Float64}
         Km_Min      ::	Vector{Float64}
         Φ_Min       ::	Vector{Float64}
			
         θs_Max      ::	Vector{Float64}
         θr_Max      ::	Vector{Float64}
         N_Max       ::	Vector{Float64}
         Ψvg_Max     ::	Vector{Float64}
         ΨmacMat_Max :: Vector{Float64}
         Ks_Max      ::	Vector{Float64}
         Km_Max      ::	Vector{Float64}
         Φ_Max       ::	Vector{Float64}
		end # struct VANGENUCHTEN


		Base.@kwdef mutable struct BROOKS_COREY # <>=<>=<>=<>=<>=<>=<>=<>=<>
         θs          ::	Vector{Float64}
         θr          ::	Vector{Float64}
         λbc         ::	Vector{Float64}
         Ψbc         ::	Vector{Float64}
         ΨmacMat     :: Vector{Float64}
         Ks          ::	Vector{Float64}
         Φ           ::	Vector{Float64}
         Ψga         ::	Vector{Float64}
			
         θs_Min      ::	Vector{Float64}
         θr_Min      ::	Vector{Float64}
         λbc_Min     ::	Vector{Float64}
         Ψbc_Min     ::	Vector{Float64}
         ΨmacMat_Min :: Vector{Float64}
         Ks_Min      ::	Vector{Float64}
         Φ_Min       ::	Vector{Float64}
			
         θs_Max      ::	Vector{Float64}
         θr_Max      ::	Vector{Float64}
         λbc_Max     ::	Vector{Float64}
         Ψbc_Max     ::	Vector{Float64}
         ΨmacMat_Max :: Vector{Float64}
         Ks_Max      ::	Vector{Float64}
         Φ_Max       ::	Vector{Float64}
		end # struct BROOKS COREY


		Base.@kwdef mutable struct CLAPP_HORNBERGER # <>=<>=<>=<>=<>=<>=<>=<>=<>
         θs          :: Vector{Float64}
         θr          :: Vector{Float64}
         λch         :: Vector{Float64}
         Ψch         :: Vector{Float64}
         ΨmacMat     :: Vector{Float64}
         Ks          :: Vector{Float64}
         Φ           :: Vector{Float64}
         Ψga         :: Vector{Float64}
			
         θs_Min      :: Vector{Float64}
         θr_Min      :: Vector{Float64}
         λch_Min     :: Vector{Float64}
         Ψch_Min     :: Vector{Float64}
         ΨmacMat_Min :: Vector{Float64}
         Ks_Min      :: Vector{Float64}
         Φ_Min       :: Vector{Float64}
			
         θs_Max      :: Vector{Float64}
         θr_Max      :: Vector{Float64}
         λch_Max     :: Vector{Float64}
         Ψch_Max     :: Vector{Float64}
         ΨmacMat_Max :: Vector{Float64}
         Ks_Max      :: Vector{Float64}
         Φ_Max       :: Vector{Float64}
		end # struct CLAPP_HORNBERGER


		Base.@kwdef mutable struct HYDRO_OTHER # <>=<>=<>=<>=<>=<>=<>=<>=<>
         Macro_Perc                         :: Vector{Float64}
         Nse                                :: Vector{Float64}
         Nse_θΨ                             :: Vector{Float64}
         Nse_KΨ                             :: Vector{Float64}
         NseWilmot_θΨ                       :: Vector{Float64}
         NseWilmot_KΨ                       :: Vector{Float64}
         Rmse                               :: Vector{Float64}
         Rmse_θΨ                            :: Vector{Float64}
         Rmse_KΨ                            :: Vector{Float64}
         NseConcordanceCorelationCoeficient :: Vector{Float64}
		end

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROSTRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDROSTRUCT(optionₘ, NiZ::Int64)
			# For all models
         θs          = fill(0.0::Float64, NiZ)
         θr          = fill(0.0::Float64, NiZ)
         Ks          = fill(0.0::Float64, NiZ)
         θsMacMat    = fill(0.0::Float64, NiZ)
         θsMacMat_ƞ  = fill(0.0::Float64, NiZ)
         ΨmacMat     = fill(100.0::Float64, NiZ)
         Φ           = fill(0.0::Float64, NiZ)
         So          = fill(0.0::Float64, NiZ)

         θs_Min      = fill(0.0::Float64, NiZ)
         θr_Min      = fill(0.0::Float64, NiZ)
         Ks_Min      = fill(0.0::Float64, NiZ)
         ΨmacMat_Min = fill(0.0::Float64, NiZ)
         Φ_Min       = fill(0.0::Float64, NiZ)
         So_Min      = fill(0.0::Float64, NiZ)
			
         θs_Max      = fill(0.0::Float64, NiZ)
         θr_Max      = fill(0.0::Float64, NiZ)
         Ks_Max      = fill(0.0::Float64, NiZ)
         ΨmacMat_Max = fill(0.0::Float64, NiZ)
         Φ_Max       = fill(0.0::Float64, NiZ)
         So_Max      = fill(0.0::Float64, NiZ)
			
			if optionₘ.HydroModel⍰ == "Kosugi" # <>=<>=<>=<>=<>
            σ              = fill(0.0::Float64, NiZ)
            Ψm             = fill(0.0::Float64, NiZ)
            σMac           = fill(0.0::Float64, NiZ)
            ΨmMac          = fill(0.0::Float64, NiZ)
            τa             = fill(0.0::Float64, NiZ)
            τb             = fill(0.0::Float64, NiZ)
            τc             = fill(0.0::Float64, NiZ)
            τₚ             = fill(0.0::Float64, NiZ)
            τaMac          = fill(0.0::Float64, NiZ)
            τbMac          = fill(0.0::Float64, NiZ)
            τcMac          = fill(0.0::Float64, NiZ)
				
            Ψm_Min         = fill(0.0::Float64, NiZ)
            σ_Min          = fill(0.0::Float64, NiZ)
            θsMacMat_ƞ_Min = fill(0.0::Float64, NiZ)
            σMac_Min       = fill(0.0::Float64, NiZ)
            ΨmMac_Min      = fill(0.0::Float64, NiZ)
            θsMacMat_Min   = fill(0.0::Float64, NiZ)
            τa_Min         = fill(0.0::Float64, NiZ)
            τb_Min         = fill(0.0::Float64, NiZ)
            τc_Min         = fill(0.0::Float64, NiZ)
            τₚ_min         = fill(0.0::Float64, NiZ)
            τaMac_Min      = fill(0.0::Float64, NiZ)
            τbMac_Min      = fill(0.0::Float64, NiZ)
            τcMac_Min      = fill(0.0::Float64, NiZ)

            Ψm_Max         = fill(0.0::Float64, NiZ)
            σ_Max          = fill(0.0::Float64, NiZ)
            θsMacMat_ƞ_Max = fill(0.0::Float64, NiZ)
            σMac_Max       = fill(0.0::Float64, NiZ)
            ΨmMac_Max      = fill(0.0::Float64, NiZ)
            θsMacMat_Max   = fill(0.0::Float64, NiZ)
            τa_Max         = fill(0.0::Float64, NiZ)
            τb_Max         = fill(0.0::Float64, NiZ)
            τc_Max         = fill(0.0::Float64, NiZ)
            τₚ_max         = fill(0.0::Float64, NiZ)
            τaMac_Max      = fill(0.0::Float64, NiZ)
            τbMac_Max      = fill(0.0::Float64, NiZ)
            τcMac_Max      = fill(0.0::Float64, NiZ)
          
				hydro = KOSUGI(θs, θr, Ks, σ, Ψm, θsMacMat_ƞ, σMac, ΨmMac, ΨmacMat, So, θsMacMat, Φ,  τa, τb, τc, τₚ, τaMac, τbMac, τcMac,θs_Min, θr_Min, Ks_Min, σ_Min,Ψm_Min, θsMacMat_ƞ_Min, σMac_Min, ΨmMac_Min, ΨmacMat_Min, So_Min, θsMacMat_Min, Φ_Min,  τa_Min, τb_Min, τc_Min, τₚ_min, τaMac_Min, τbMac_Min, τcMac_Min, θs_Max, θr_Max, Ks_Max, σ_Max, Ψm_Max, θsMacMat_ƞ_Max, σMac_Max, ΨmMac_Max, ΨmacMat_Max, So_Max, θsMacMat_Max, Φ_Max, τa_Max, τb_Max, τc_Max, τₚ_max, τaMac_Max, τbMac_Max, τcMac_Max)   
           
				return hydro

			elseif optionₘ.HydroModel⍰=="Vangenuchten" || optionₘ.HydroModel⍰=="VangenuchtenJules" # <>=<>=<>=<>=<>
            N       = fill(0.0::Float64, NiZ)
            Ψvg     = fill(0.0::Float64, NiZ)
            Km      = fill(0.0::Float64, NiZ)

            θs_Min  = fill(0.0::Float64, NiZ)
            θr_Min  = fill(0.0::Float64, NiZ)
            N_Min   = fill(0.0::Float64, NiZ)
            Ψvg_Min = fill(0.0::Float64, NiZ)
            Ks_Min  = fill(0.0::Float64, NiZ)
            Km_Min  = fill(0.0::Float64, NiZ)
			
            θs_Max  = fill(0.0::Float64, NiZ)
            θr_Max  = fill(0.0::Float64, NiZ)
            N_Max   = fill(0.0::Float64, NiZ)
            Ψvg_Max = fill(0.0::Float64, NiZ)
            Ks_Max  = fill(0.0::Float64, NiZ)
            Km_Max  = fill(0.0::Float64, NiZ)

				hydro = VANGENUCHTEN(θs, θr, N, Ψvg, ΨmacMat, Ks, Km, Φ, θs_Min,θr_Min, N_Min, Ψvg_Min, ΨmacMat_Min, Ks_Min, Km_Min, Φ_Min, θs_Max, θr_Max, N_Max, Ψvg_Max, ΨmacMat_Max, Ks_Max, Km_Max, Φ_Max) 
				return hydro


			elseif optionₘ.HydroModel⍰ == "BrooksCorey" # <>=<>=<>=<>=<>=<>
            λbc     = fill(0.0::Float64, NiZ)
            Ψbc     = fill(0.0::Float64, NiZ)
            λbc_Min = fill(0.0::Float64, NiZ)
            Ψbc_Min = fill(0.0::Float64, NiZ)
            λbc_Max = fill(0.0::Float64, NiZ)
            Ψbc_Max = fill(0.0::Float64, NiZ)
            Ψga     = fill(0.0::Float64, NiZ)

				hydro = BROOKS_COREY(θs, θr, λbc, Ψbc, ΨmacMat, Ks, Φ, Ψga, θs_Min, θr_Min, λbc_Min, Ψbc_Min, ΨmacMat_Min, Ks_Min, Φ_Min, θs_Max, θr_Max, λbc_Max, Ψbc_Max, ΨmacMat_Max, Ks_Max, Φ_Max)
				return hydro

			elseif optionₘ.HydroModel⍰ == "ClappHornberger" # <>=<>=<>=<>=<>=<>
				λch = fill(0.0::Float64, NiZ)
				Ψch  = fill(0.0::Float64, NiZ)
            Ψga  = fill(0.0::Float64, NiZ)

				λch_Min = fill(0.0::Float64, NiZ)
				Ψch_Min  = fill(0.0::Float64, NiZ)

				λch_Max = fill(0.0::Float64, NiZ)
				Ψch_Max  = fill(0.0::Float64, NiZ)

				hydro = CLAPP_HORNBERGER(θs, θr, λch, Ψch, ΨmacMat, Ks, Φ, Ψga, θs_Min, θr_Min, λch_Min, Ψch_Min, ΨmacMat_Min, Ks_Min, Φ_Max, θs_Max, θr_Max, λch_Max, Ψch_Max, ΨmacMat_Max, Ks_Max, Φ_Max)
				return hydro
			end # optionₘ.HydroModel⍰

		end #  function HYDROSTRUCT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYDRO_OTHER
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYDRO_OTHERS(NiZ::Int64)
            Macro_Perc                         = fill(0.0::Float64, NiZ)
            Nse                                = fill(0.0::Float64, NiZ)
            Nse_θΨ                             = fill(0.0::Float64, NiZ)
            Nse_KΨ                             = fill(0.0::Float64, NiZ)
            NseWilmot_θΨ                       = fill(0.0::Float64, NiZ)
            NseWilmot_KΨ                       = fill(0.0::Float64, NiZ)
            Rmse                               = fill(0.0::Float64, NiZ)
            Rmse_θΨ                            = fill(0.0::Float64, NiZ)
            Rmse_KΨ                            = fill(0.0::Float64, NiZ)
            NseConcordanceCorelationCoeficient = fill(0.0::Float64, NiZ)

				hydroOther = HYDRO_OTHER(Macro_Perc, Nse, Nse_θΨ, Nse_KΨ, NseWilmot_θΨ, NseWilmot_KΨ, Rmse, Rmse_θΨ, Rmse_KΨ, NseConcordanceCorelationCoeficient)	
				return hydroOther
			end  # function: HYDRO_OTHER
end  # module: hydroStruct

# ............................................................