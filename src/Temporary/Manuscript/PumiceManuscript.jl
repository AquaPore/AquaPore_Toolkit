
Path =  "D:\\MAIN\\MODELS\\AquaPore_Toolkit\\src\\"

# include(Path * "Including.jl")
include(Path * "Cst.jl")
include(Path * "hydro//HydroRelation.jl")
include(Path * "hydro//Wrc.jl")
include(Path * "hydro//Kunsat.jl")

module pumiceManuscript
	using CairoMakie, ColorSchemes
	import SpecialFunctions: erfc, erfcinv

	import ..wrc, ..kunsat, ..hydroRelation, ..cst

	export PLOTTING_PORESIZE

	Pσ₁ = 1.0
	Pσ₂ = 2.0
	Pσ₃ = 3.0
	Pσ₄ = 4.0


	# ------------------------------------------------------------------

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HydroModels
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDRO_MODELS(;θs, θsMacMat_η, θr=0.0, σ, Ks, ΨmacMat, τa=0.5, τb, τc, τaMac=0.5, τbMac, τcMac, KosugiModel_θΨ⍰, Option_KosugiModel_KΨ⍰, Pσ=3, Pσ_Mac=3)

			θsMacMat   = (θs - θr) * θsMacMat_η + θr

			# Deriving macropore hydraulic parameters from ΨmacMat
				σMac    = hydroRelation.FUNC_ΨmacMat_2_σMac(;ΨmacMat=ΨmacMat, Pσ_Mac=3)
				ΨmMac   = hydroRelation.FUNC_ΨmacMat_2_ΨmMac(;ΨmacMat=ΨmacMat, σMac)

				Ψm_Min  = hydroRelation.FUNC_σ_2_Ψm(;ΨmacMat= ΨmacMat, σ, Pσ=Pσ₃, 🎏_Min=true)
				Ψm_Max  = hydroRelation.FUNC_σ_2_Ψm(;ΨmacMat=ΨmacMat, σ, Pσ=Pσ₃, 🎏_Min=false)

				Tb_Max = 2.0; Tc_Max = 4.0
				Tb    = Tb_Max * (1.0 - τb)
				TbMac = Tb_Max * (1.0 - τbMac)
				Tc    = Tc_Max * (1.0 - τc)
				# Tc= σ ^ -0.59
				TcMac = Tc_Max * (1.0 - τcMac)
			
				KsMac_Min, KsMat_Min = kunsat.kg.KS_MATMAC_ΨmacMat(θs, θsMacMat, θr, Ψm_Min, σ, ΨmMac, σMac, Ks, Tb, Tc, TbMac, TcMac, KosugiModel_θΨ⍰)
				KsMac_Max, KsMat_Max = kunsat.kg.KS_MATMAC_ΨmacMat(θs, θsMacMat, θr, Ψm_Max, σ, ΨmMac, σMac, Ks, Tb, Tc, TbMac, TcMac, KosugiModel_θΨ⍰)

			#  For every ψ
            Ψ_Min_Log = log10(0.0001)
            Ψ_Max_Log = log10(1500_00)

				Ψ = 10.0.^(collect(Ψ_Min_Log:0.0001:Ψ_Max_Log))
				N = length(Ψ)

			# Initiating
            Kunsat_Min = Array{Float64}(undef, N)
            Kunsat_Max = Array{Float64}(undef, N)
            θDual_Min  = Array{Float64}(undef, N)
            θDual_Max  = Array{Float64}(undef, N)

			# Looping
			for iΨ=1:N
		
				# ∂θ∂Ψ_1[iΨ] = wrc.kg.∂θ∂Ψ_NORM(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Min, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac)

				# ∂θ∂Ψ_2[iΨ] =  wrc.kg.∂θ∂Ψ_NORM(Ψ₁=Ψ[iΨ], θs=θs, θr=θr,  Ψm=Ψm_Max, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac)

				θDual_Min[iΨ] = wrc.kg.Ψ_2_θ(;Ψ₁=Ψ[iΨ], θs, θsMacMat, θr, Ψm=Ψm_Min, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰)

				
				θDual_Max[iΨ] = wrc.kg.Ψ_2_θ(;Ψ₁=Ψ[iΨ], θs, θsMacMat, θr, Ψm=Ψm_Max, σ, ΨmMac, ΨmacMat, σMac, KosugiModel_θΨ⍰)

				Kunsat_Min[iΨ] = kunsat.kg.KUNSAT_θΨSe(;Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Min, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac, Ks=Ks, τa=τa, τb=τb, τc=τc, τaMac=τaMac, τbMac=τbMac, τcMac=τcMac, Option_KosugiModel_KΨ⍰, KosugiModel_θΨ⍰)

				Kunsat_Max[iΨ] = kunsat.kg.KUNSAT_θΨSe(;Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Max, σ=σ, θsMacMat=θsMacMat,ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac, Ks=Ks, τa=τa, τb=τb, τc=τc, τaMac=τaMac, τbMac=τbMac, τcMac=τcMac, Option_KosugiModel_KΨ⍰, KosugiModel_θΨ⍰)
			end
			
		return Ks, KsMac_Max, KsMac_Max, KsMac_Min, KsMac_Min, KsMat_Max, KsMat_Min, Kunsat_Max, Kunsat_Min, N, Pσ_Mac, θDual_Max, θDual_Min, θs, θsMacMat, σ, Ψ, ΨmacMat
		end  # function: HydroModels
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOTTING_PORESIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PLOTTING_PORESIZE()

		# Plotting parameters
         ColourOption_No    = 8
         Linewidth          = 2
         height             = 200
         labelsize          = 12
         textcolor          = :blue
         textsize           = 17
         titlecolor         = :darkgoldenrod4
         titlesize          = 20.0
         width              = height * 3.0
         xgridstyle         = :dash
         xgridvisible       = true
         xlabelSize         = 20
         xlabelpadding      = 5
         xminortickalign    = 1.0
         xminorticksvisible = true
         xtickalign         = 0.9 # 0 is inside and 1 is outside
         xticklabelrotation = π/4
         xticksize          = 10
         xticksmirrored     = false
         xtickwidt          = 0.5
         xtrimspine         = false
         ygridstyle         = :dash
         ygridvisible       = false
         ylabelpadding      = xlabelpadding
         ylabelsize         = xlabelSize
         yminortickalign    = xminortickalign
         yminorticksvisible = true
         ytickalign         = xtickalign
         yticksize          = xticksize
         yticksmirrored     = false
         ytickwidt          = xtickwidt
         ytrimspine         = false

			ColourOption = [:glasbey_hv_n256,:seaborn_bright,:seaborn_colorblind,:seaborn_dark,:seaborn_deep,:tab10,:tableau_10,:tol_bright]

			Colormap = cgrad(colorschemes[ColourOption[ColourOption_No]], size(colorschemes[ColourOption[ColourOption_No]]), categorical = true)

			Ψticks = [0, 50, 100, 500, 1000,5000,100_00, 500_00, 1000_00, 2000_00] # mm

		# Starting to plot	
			CairoMakie.activate!(type="svg", pt_per_unit=1)
			Fig =  Figure(figure_padding = 10; fonts = ( ; regular="CMU Serif")) 


		Axis_θΨ = Axis(Fig[1,1], xlabel= L"$ψ$ [kPa]", ylabel=L"$θ$ [L³ L⁻³]", title="Sandy soil", titlecolor=titlecolor, xticklabelrotation=xticklabelrotation, ylabelsize=ylabelsize, xlabelsize=xlabelSize, xticksize=xticksize, yticksize=yticksize, width=width, height=height,   titlesize=titlesize,  xgridvisible=xgridvisible, ygridvisible=ygridvisible, xminorticksvisible=xminorticksvisible, yminorticksvisible=yminorticksvisible, xtickwidth=xtickwidt, ytickwidth=ytickwidt, xtickalign=xtickalign, ytickalign=ytickalign, xticksmirrored=xticksmirrored, yticksmirrored=yticksmirrored,  xtrimspine=xtrimspine,  ytrimspine=ytrimspine, xgridstyle=xgridstyle, ygridstyle=ygridstyle, yminorticks=IntervalsBetween(5), xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding, xminortickalign=xminortickalign, yminortickalign=yminortickalign)

		Label(Fig[1, 1, TopRight()], "(A1)", fontsize=18, font=:bold, padding=(-50, 5, -50, 5), halign=:right)

		Ks, KsMac_Max, KsMac_Max, KsMac_Min, KsMac_Min, KsMat_Max, KsMat_Min, Kunsat_Max, Kunsat_Min, N, Pσ_Mac, θDual_Max, θDual_Min, θs, θsMacMat, σ, Ψ, ΨmacMat =HYDRO_MODELS(θs=0.5, θsMacMat_η=0.75, θr=0.0, σ=1.0, Ks=0.08, ΨmacMat=100.0, τa=0.5, τb=0.6, τc=0.5, τaMac=0.5, τbMac=0.9, τcMac=0.03, KosugiModel_θΨ⍰="ΨmacMat", Option_KosugiModel_KΨ⍰="ΨmacMat")

				Label(Fig[1, 1:2, Top()], "Functions θ(ψ) & K(ψ) macropore Pσ_Mac=$Pσ_Mac", valign=:bottom, font=:bold, padding=(0, 0,50, 0), color=:navajowhite4,  fontsize=titlesize*1.2)
					
				Ψ_Log = Array{Float64}(undef, N)
				for iZ=1:N
					Ψ_Log[iZ] = log1p(Ψ[iZ])
				end

				Axis_θΨ.xticks = (log1p.(Ψticks), string.(cst.Mm_2_kPa .* Ψticks))
				hidexdecorations!(Axis_θΨ, ticks=false, grid=false)

				lines!(Axis_θΨ, Ψ_Log, θDual_Min, linewidth=Linewidth, label="Macro σ =$σ, ψₘ_Min", color=Colormap[1], linestyle=:dashdot)
				lines!(Axis_θΨ, Ψ_Log, θDual_Max, linewidth=Linewidth, label="Macro σ =$σ, ψₘ_Max",  color=Colormap[1])

				lines!(Axis_θΨ,[Point(log1p(ΨmacMat),0), Point(log1p(ΨmacMat), θs)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
				lines!(Axis_θΨ,[Point(log1p(0), θsMacMat ), Point(log1p(ΨmacMat), θsMacMat)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
				lines!(Axis_θΨ,[Point(log1p(0), θs), Point(log1p(ΨmacMat), θs)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)

				text!(log1p(0), θsMacMat, text =L"θ_{sMacMat}", align=(:left,:top),  color=textcolor, fontsize=textsize)
				text!(log1p(0), θs, text =L"θ_{s}", align=(:left,:top),  color=textcolor, fontsize=textsize)
				text!(log1p(ΨmacMat), 0, text =L"ψ_{macMat}", align=(:left,:bottom), rotation = π/2,  color=textcolor, fontsize=textsize)
				
				Ks_Trad, KsMac_Max_Trad, KsMac_Max_trad, KsMac_Min_Trad, KsMac_Min_Trad, KsMat_Max_Trad, KsMat_Min_trad, Kunsat_Max_Trad, Kunsat_Min_Trad, N, Pσ_Mac, θDual_Max_Trad, θDual_Min_Trad, θs, θsMacMat, σ, Ψ, ΨmacMat = HYDRO_MODELS(;θs=0.5, θsMacMat_η=0.75, θr=0.0, σ=1.0, Ks=0.08, ΨmacMat=100.0, τa=0.5, τb=0.6, τc=0.6, τaMac=0.5, τbMac=0.6, τcMac=0.0, KosugiModel_θΨ⍰="Traditional", Option_KosugiModel_KΨ⍰="Traditional")

				lines!(Axis_θΨ, Ψ_Log, θDual_Min_Trad, linewidth=Linewidth, label="Trad σ =$σ, ψₘ_Min", color=Colormap[2], linestyle=:dashdot)
				lines!(Axis_θΨ, Ψ_Log, θDual_Max_Trad, linewidth=Linewidth, label="Trad σ =$σ, ψₘ_Max",  color=Colormap[2])


		Axis_KΨ = Axis(Fig[2,1], xlabel=L"$ψ$ [kPa]", ylabel=L"$K(ψ)$ [cm h⁻¹]", xticklabelrotation=xticklabelrotation, ylabelsize=ylabelsize,	xlabelsize=xlabelSize, xticksize=xticksize, yticksize=yticksize, width=width, height=height, xgridvisible=xgridvisible, ygridvisible=ygridvisible, xminorticksvisible=xminorticksvisible, yminorticksvisible=yminorticksvisible, xtickwidth=xtickwidt, ytickwidth=ytickwidt, xticksmirrored=xticksmirrored,  yticksmirrored=yticksmirrored, xtrimspine=xtrimspine,  ytrimspine=ytrimspine, xtickalign=xtickalign, ytickalign=ytickalign, yminorticks=IntervalsBetween(5), xgridstyle=xgridstyle, ygridstyle=ygridstyle, xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding, xminortickalign=xminortickalign, yminortickalign=yminortickalign)

			Label(Fig[2, 1, TopRight()], "(A2)", fontsize=18, font=:bold, padding = (-50, 5, -50, 5), halign=:right)

			Axis_KΨ.xticks = (log1p.(Ψticks), string.(cst.Mm_2_kPa .* Ψticks))

			lines!(Axis_KΨ, Ψ_Log , cst.MmS_2_CmH .* Kunsat_Min,  linewidth=Linewidth, color=Colormap[1], linestyle=:dashdot)
			lines!(Axis_KΨ, Ψ_Log , cst.MmS_2_CmH .* Kunsat_Max, linewidth=Linewidth, color=Colormap[1])

			lines!(Axis_KΨ, Ψ_Log , cst.MmS_2_CmH .* Kunsat_Min_Trad,  linewidth=Linewidth, color=Colormap[2], linestyle=:dashdot)
			lines!(Axis_KΨ, Ψ_Log , cst.MmS_2_CmH .* Kunsat_Max_Trad, linewidth=Linewidth, color=Colormap[2])

			lines!(Axis_KΨ,[Point(log1p(ΨmacMat),0), Point(log1p(ΨmacMat), cst.MmS_2_CmH .* Ks)], color=:navyblue, linewidth=Linewidth/2.0,  linestyle=:dash)
			lines!(Axis_KΨ,[Point(log1p(0), cst.MmS_2_CmH .* KsMat_Min ), Point(log1p(ΨmacMat), cst.MmS_2_CmH .*  KsMat_Min)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
			lines!(Axis_KΨ,[Point(log1p(0), cst.MmS_2_CmH .* KsMat_Max ), Point(log1p(ΨmacMat), cst.MmS_2_CmH .*  KsMat_Max)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
			lines!(Axis_KΨ,[Point(log1p(0),  cst.MmS_2_CmH .* Ks), Point(log1p(ΨmacMat), cst.MmS_2_CmH .* Ks)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)

			text!(log1p(0),  cst.MmS_2_CmH .* KsMat_Min, text =L"K_{sMacMin}", align=(:left,:top), color=textcolor, fontsize=textsize)
			text!(log1p(0),  cst.MmS_2_CmH .* KsMat_Max, text =L"K_{sMacMax}", align=(:left,:top), color=textcolor, fontsize=textsize)
			text!(log1p(0),  cst.MmS_2_CmH .* Ks, text =L"K_{s}", align=(:left,:top), color=textcolor, fontsize=textsize)
			text!(log1p(ΨmacMat), 0, text =L"ψ_{macMat}", align=(:left,:bottom), rotation = π/2, color=textcolor, fontsize=textsize)

			# Legend
				Legend(Fig[3,1], Axis_θΨ, framecolor=(:grey, 0.5), labelsize=labelsize, valign=:top, padding=5, tellheight=true, tellwidt=true, nbanks=2, backgroundcolor = (:grey90, 0.25))

		# ========================================

		Axis_θΨ2 = Axis(Fig[1,2], xlabel= L"$ψ$ [kPa]", ylabel=L"$θ$ [L³ L⁻³]", title="Clay soils",  titlecolor=titlecolor, xticklabelrotation=xticklabelrotation, ylabelsize=ylabelsize, xlabelsize=xlabelSize, xticksize=xticksize, yticksize=yticksize, width=width, height=height, titlesize=titlesize,  xgridvisible=xgridvisible, ygridvisible=ygridvisible, xminorticksvisible=xminorticksvisible, yminorticksvisible=yminorticksvisible, xtickwidth=xtickwidt, ytickwidth=ytickwidt, xtickalign=xtickalign, ytickalign=ytickalign, xticksmirrored=xticksmirrored, yticksmirrored=yticksmirrored,  xtrimspine=xtrimspine,  ytrimspine=ytrimspine, xgridstyle=xgridstyle, ygridstyle=ygridstyle, yminorticks=IntervalsBetween(5), xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding, xminortickalign=xminortickalign, yminortickalign=yminortickalign)

			Label(Fig[1, 2, TopRight()], "(B1)", fontsize=18, font=:bold, padding=(-50, 5, -50, 5), halign=:right)

			Ks, KsMac_Max, KsMac_Max, KsMac_Min, KsMac_Min, KsMat_Max, KsMat_Min, Kunsat_Max, Kunsat_Min, N, Pσ_Mac, θDual_Max, θDual_Min, θs, θsMacMat, σ, Ψ, ΨmacMat =HYDRO_MODELS(θs=0.5, θsMacMat_η=0.75, θr=0.0, σ=3.0, Ks=0.08, ΨmacMat=100.0, τa=0.5, τb=0.5, τc=0.6, τaMac=0.5, τbMac=0.6, τcMac=0.0, KosugiModel_θΨ⍰="ΨmacMat", Option_KosugiModel_KΨ⍰="ΨmacMat")

				Ψ_Log = Array{Float64}(undef, N)
				for iZ=1:N
					Ψ_Log[iZ] = log1p(Ψ[iZ])
				end

				Axis_θΨ2.xticks = (log1p.(Ψticks), string.(cst.Mm_2_kPa .* Ψticks))
				hidexdecorations!(Axis_θΨ2, ticks=false, grid=false)

				lines!(Axis_θΨ2, Ψ_Log, θDual_Min, linewidth=Linewidth, label="Macro σ =$σ, ψₘ_Min", color=Colormap[3], linestyle=:dashdot)
				lines!(Axis_θΨ2, Ψ_Log, θDual_Max, linewidth=Linewidth, label="Macro σ =$σ, ψₘ_Max",  color=Colormap[3])

				lines!(Axis_θΨ2,[Point(log1p(ΨmacMat),0), Point(log1p(ΨmacMat), θs)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
				lines!(Axis_θΨ2,[Point(log1p(0), θsMacMat ), Point(log1p(ΨmacMat), θsMacMat)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
				lines!(Axis_θΨ2,[Point(log1p(0), θs), Point(log1p(ΨmacMat), θs)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)

				text!(log1p(0), θsMacMat, text =L"θ_{sMacMat}", align=(:left,:top),  color=textcolor, fontsize=textsize)
				text!(log1p(0), θs, text =L"θ_{s}", align=(:left,:top),  color=textcolor, fontsize=textsize)
				text!(log1p(ΨmacMat), 0, text =L"ψ_{macMat}", align=(:left,:bottom), rotation = π/2,  color=textcolor, fontsize=textsize)
				
				Ks_Trad, KsMac_Max_Trad, KsMac_Max_trad, KsMac_Min_Trad, KsMac_Min_Trad, KsMat_Max_Trad, KsMat_Min_trad, Kunsat_Max_Trad, Kunsat_Min_Trad, N, Pσ_Mac, θDual_Max_Trad, θDual_Min_Trad, θs, θsMacMat, σ, Ψ, ΨmacMat = HYDRO_MODELS(;θs=0.5, θsMacMat_η=0.75, θr=0.0, σ=3.0, Ks=0.08, ΨmacMat=100.0, τa=0.5, τb=0.6, τc=0.6, τaMac=0.5, τbMac=0.6, τcMac=0.0, KosugiModel_θΨ⍰="Traditional", Option_KosugiModel_KΨ⍰="Traditional")

				lines!(Axis_θΨ2, Ψ_Log, θDual_Min_Trad, linewidth=Linewidth, label="Trad σ =$σ, ψₘ_Min", color=Colormap[4], linestyle=:dashdot)
				lines!(Axis_θΨ2, Ψ_Log, θDual_Max_Trad, linewidth=Linewidth, label="Trad σ =$σ, ψₘ_Max",  color=Colormap[4])


		Axis_KΨ2 = Axis(Fig[2,2], xlabel=L"$ψ$ [kPa]", ylabel=L"$K(ψ)$ [cm h⁻¹]", xticklabelrotation=xticklabelrotation, ylabelsize=ylabelsize,	xlabelsize=xlabelSize, xticksize=xticksize, yticksize=yticksize, width=width, height=height, xgridvisible=xgridvisible, ygridvisible=ygridvisible, xminorticksvisible=xminorticksvisible, yminorticksvisible=yminorticksvisible, xtickwidth=xtickwidt, ytickwidth=ytickwidt, xticksmirrored=xticksmirrored,  yticksmirrored=yticksmirrored, xtrimspine=xtrimspine,  ytrimspine=ytrimspine, xtickalign=xtickalign, ytickalign=ytickalign, yminorticks=IntervalsBetween(5), xgridstyle=xgridstyle, ygridstyle=ygridstyle, xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding, xminortickalign=xminortickalign, yminortickalign=yminortickalign )

			Label(Fig[2, 2, TopRight()], "(B2)", fontsize=18, font=:bold, padding = (-50, 5, -50, 5), halign=:right)

			Axis_KΨ2.xticks = (log1p.(Ψticks), string.(cst.Mm_2_kPa .* Ψticks))

			lines!(Axis_KΨ2, Ψ_Log , cst.MmS_2_CmH .* Kunsat_Min,  linewidth=Linewidth, color=Colormap[3], linestyle=:dashdot)
			lines!(Axis_KΨ2, Ψ_Log , cst.MmS_2_CmH .* Kunsat_Max, linewidth=Linewidth, color=Colormap[3])

			lines!(Axis_KΨ2, Ψ_Log , cst.MmS_2_CmH .* Kunsat_Min_Trad,  linewidth=Linewidth, color=Colormap[4], linestyle=:dashdot)
			lines!(Axis_KΨ2, Ψ_Log , cst.MmS_2_CmH .* Kunsat_Max_Trad, linewidth=Linewidth, color=Colormap[4])

			lines!(Axis_KΨ2,[Point(log1p(ΨmacMat),0), Point(log1p(ΨmacMat), cst.MmS_2_CmH .* Ks)], color=:navyblue, linewidth=Linewidth/2.0,  linestyle=:dash)
			lines!(Axis_KΨ2,[Point(log1p(0), cst.MmS_2_CmH .* KsMat_Min ), Point(log1p(ΨmacMat), cst.MmS_2_CmH .*  KsMat_Min)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
			lines!(Axis_KΨ2,[Point(log1p(0), cst.MmS_2_CmH .* KsMat_Max ), Point(log1p(ΨmacMat), cst.MmS_2_CmH .*  KsMat_Max)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
			lines!(Axis_KΨ2,[Point(log1p(0),  cst.MmS_2_CmH .* Ks), Point(log1p(ΨmacMat), cst.MmS_2_CmH .* Ks)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)

			text!(log1p(0),  cst.MmS_2_CmH .* KsMat_Min, text =L"K_{sMacMin}", align=(:left,:top), color=textcolor, fontsize=textsize)
			text!(log1p(0),  cst.MmS_2_CmH .* KsMat_Max, text =L"K_{sMacMax}", align=(:left,:top), color=textcolor, fontsize=textsize)
			text!(log1p(0),  cst.MmS_2_CmH .* Ks, text =L"K_{s}", align=(:left,:top), color=textcolor, fontsize=textsize)
			text!(log1p(ΨmacMat), 0, text =L"ψ_{macMat}", align=(:left,:bottom), rotation = π/2, color=textcolor, fontsize=textsize)

			# Legend
				Legend(Fig[3,2], Axis_θΨ2, framecolor=(:grey, 0.5), labelsize=labelsize, valign=:top, padding=5, tellheight=true, tellwidt=true, nbanks=2, backgroundcolor = (:grey90, 0.25))

			# General
				resize_to_layout!(Fig)
				trim!(Fig.layout)
				colgap!(Fig.layout, 10)
				rowgap!(Fig.layout, 10)

				Path = raw"D:\TEMP\Plots\MacroporeThetaH.svg"
				save(Path, Fig)
				display(Fig)

		
	return nothing
	end  # function: PLOTTING_PORESIZE



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOTTING_KUNSAT_MACRO
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PLOTTING_KUNSAT_MACRO()

		function RELATIONSHIPS_MAC(ΨmacMat; Pσ_Mac=3)
			σMac    = hydroRelation.FUNC_ΨmacMat_2_σMac(;ΨmacMat, Pσ_Mac)
			ΨmMac   = hydroRelation.FUNC_ΨmacMat_2_ΨmMac(;ΨmacMat=ΨmacMat, σMac)
		return σMac, ΨmMac
		end  # function: RELATIONSHIPS_MAC

		function RELATIONSHIPS_MAT(ΨmacMat, σ; Pσ=3)
			ΨmacMat₂ = exp((log(√ΨmacMat) + log(ΨmacMat)) * 0.5)
			Ψm  = hydroRelation.FUNC_σ_2_Ψm(;ΨmacMat=ΨmacMat₂, σ, Pσ=Pσ₃, 🎏_Min=false)
		return Ψm
		end

		function KUNSAT_MAT_η(;Ψ₁=Ψ₁, θs=1.0, θsMacMat=0.8, θr=0.0, σ, ΨmacMat, τb, τbMac, Ks=1.0, τa=0.5, τaMac=0.5, τc=1.0, τcMac=2.0)
			Ψm = RELATIONSHIPS_MAT(ΨmacMat, σ)
			σMac, ΨmMac = RELATIONSHIPS_MAC(ΨmacMat)

			Kunsat_Mat_Norm = kunsat.kg.KUNSAT_θΨSe(;Ψ₁=Ψ₁, θs, θsMacMat, θr, Ψm, σ, ΨmMac, ΨmacMat, σMac, Ks, τa, τb, τc, τaMac, τbMac, τcMac, Option_KosugiModel_KΨ⍰="ΨmacMat", KosugiModel_θΨ⍰="ΨmacMat")

			Ta, Tb, Tc, TaMac, TbMac, TcMac = kunsat.kg.TORTUOSITY(; σ, τa, τaMac, τb, τbMac, τc, τcMac)
	
			KsMac, KsMat = kunsat.kg.KS_MATMAC_ΨmacMat(θs::Float64, θsMacMat::Float64, θr::Float64, Ψm::Float64, σ::Float64, ΨmMac::Float64, σMac::Float64, Ks::Float64, Tb::Float64, Tc::Float64, TbMac::Float64, TcMac::Float64,  "ΨmacMat")

		return KsMat, Kunsat_Mat_Norm
		end

		# Parameters
         θs = 1.0
         θr = 0.0
         Ks = 1.0

		#  For every ψ
			Ψ_Min_Log = log10(0.0001)
			Ψ_Max_Log = log10(1500_00.0)
			Ψ = 10.0.^(collect(Ψ_Min_Log:0.0001:Ψ_Max_Log))
			N_Ψ  = length(Ψ)

			σ =  collect(range(0.75, stop=3.5, length=3))
			N_σ = length(σ)

			Texture = ["Sandy soils", "Silty soils", "Clay soils"] 

			Tb2 = collect(range(0.0, stop=1.0, length=5))
			N_Tb = length(Tb2)

			ΨmacMat = collect(range(60.0, stop=100.0, length=3))
			N_ΨmacMat = length(ΨmacMat)

			# θsMacMat_η = collect(range(0.75, stop=1.0, length=4))
			# N_θsMacMat_η  = length(θsMacMat_η)
      	# θsMacMat = θs .* θsMacMat_η

			# KsMac = zeros(100)


		
		# FUNCTION Kunsat_Tb
         KunsatMat_Tb = zeros(N_Ψ, N_σ, N_Tb, N_ΨmacMat)
         KsMatrice    = zeros(N_σ, N_Tb, N_ΨmacMat)
			
			for (iiΨ, iΨ) in enumerate(Ψ), (iiσ, iσ) in enumerate(σ), (iiTb, iTb) in enumerate(Tb2), (iiΨmacMat, iΨmacMat) in enumerate(ΨmacMat)

			 	KsMatrice[iiσ, iiTb, iiΨmacMat], KunsatMat_Tb[iiΨ, iiσ, iiTb, iiΨmacMat] = KUNSAT_MAT_η(;Ψ₁=iΨ, σ=iσ, τb=iTb, τbMac=iTb, ΨmacMat=iΨmacMat)
			end

		# ================================================================
				# Plotting parameters
         ColourOption_No    = 1
         Linewidth          = 2
         height             = 200
         labelsize          = 16
         textcolor          = :blue
         textsize           = 17
         titlecolor         = :navyblue
         titlesize          = 18.0
         width              = height * 3.0
         xgridstyle         = :dash
         xgridvisible       = true
         xlabelSize         = 20
         xlabelpadding      = 5
         xminortickalign    = 1.0
         xminorticksvisible = true
         xtickalign         = 0.9 # 0 is inside and 1 is outside
         xticklabelrotation = π/4
         xticksize          = 10
         xticksmirrored     = false
         xtickwidt          = 0.5
         xtrimspine         = false
         ygridstyle         = :dash
         ygridvisible       = false
         ylabelpadding      = xlabelpadding
         ylabelsize         = xlabelSize
         yminortickalign    = xminortickalign
         yminorticksvisible = true
         ytickalign         = xtickalign
         yticksize          = xticksize
         yticksmirrored     = false
         ytickwidt          = xtickwidt
         ytrimspine         = false

			ColourOption = [:glasbey_hv_n256,:seaborn_bright,:seaborn_colorblind,:seaborn_dark,:seaborn_deep,:tab10,:tableau_10,:tol_bright]

			Colormap = cgrad(colorschemes[ColourOption[ColourOption_No]], size(colorschemes[ColourOption[ColourOption_No]]), categorical = true)

			Ψticks = [0, 50, 100, 500, 1000,5000,100_00, 500_00, 1000_00, 1500_00] # mm

			Ψ_Log = Array{Float64}(undef, N_Ψ)
				for iZ=1:N_Ψ
					Ψ_Log[iZ] = log1p(Ψ[iZ])
				end

		# Starting to plot	
			CairoMakie.activate!(type="svg", pt_per_unit=1)
			Fig =  Figure(figure_padding = 10; fonts = ( ; regular="CMU Serif"), backgroundcolor = :seashell2) 

			Label(Fig[1, 1:N_σ, Top()], "Macropore parsimenous K(θ) model", valign=:bottom, font=:bold, padding=(0, 0, 50, 0), color=:darkblue,  fontsize=titlesize*1.2)

			Axis_KunsatMat_Tb = []
			for (iiσ, iσ) in enumerate(σ)
				for(iiΨmacMat, iΨmacMat) in enumerate(ΨmacMat)

				if iiσ==1
					Label(Fig[iiσ, iiΨmacMat, TopRight()], "($iiσ-$iiΨmacMat)", fontsize=18, font=:bold, padding=(-50, 5, -100, 10), halign=:right)
					else
						Label(Fig[iiσ, iiΨmacMat, TopRight()], "($iiσ-$iiΨmacMat)", fontsize=18, font=:bold, padding=(-50, 5, -50, 10), halign=:right)
					end

					Axis_KunsatMat_Tb = Axis(Fig[iiσ, iiΨmacMat], xlabel= L"$ψ$ [kPa]", ylabel=L"$K_{unsat}$ [L³ L⁻³]", title="$(Texture[iiσ]) ΨmacMat=$(Int32(floor(iΨmacMat, digits=0))) mm" ,  titlecolor=titlecolor, xticklabelrotation=xticklabelrotation, ylabelsize=ylabelsize, xlabelsize=xlabelSize, xticksize=xticksize, yticksize=yticksize, width=width, height=height, titlesize=titlesize,  xgridvisible=xgridvisible, ygridvisible=ygridvisible, xminorticksvisible=xminorticksvisible, yminorticksvisible=yminorticksvisible, xtickwidth=xtickwidt, ytickwidth=ytickwidt, xtickalign=xtickalign, ytickalign=ytickalign, xticksmirrored=xticksmirrored, yticksmirrored=yticksmirrored, xtrimspine=xtrimspine,  ytrimspine=ytrimspine, xgridstyle=xgridstyle, ygridstyle=ygridstyle, yminorticks=IntervalsBetween(5), xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding, xminortickalign=xminortickalign, yminortickalign=yminortickalign)

						Axis_KunsatMat_Tb.xticks = (log1p.(Ψticks), string.(cst.Mm_2_kPa .* Ψticks))

						if iiσ < N_σ
							hidexdecorations!(Axis_KunsatMat_Tb, ticks=false, grid=false)
						end

						if iiΨmacMat > 1
							hideydecorations!(Axis_KunsatMat_Tb, ticks=false, grid=false)
						end

						for (iiTb, iTb) in enumerate(Tb2)
							lines!(Fig[iiσ, iiΨmacMat], Ψ_Log, KunsatMat_Tb[:, iiσ, iiTb, iiΨmacMat], linewidth=Linewidth, color=Colormap[iiTb], label="σ =$(floor(σ[iiσ], digits=2)) Tb=$(floor(iTb, digits=2))")

							lines!(Axis_KunsatMat_Tb,[Point(log1p(0.0), KsMatrice[iiσ, iiTb, iiΨmacMat] ), Point(log1p(ΨmacMat[iiΨmacMat]), KsMatrice[iiσ, iiTb, iiΨmacMat])], color=Colormap[iiTb], linewidth=Linewidth/2.0, linestyle=:dash)
						end

						lines!(Axis_KunsatMat_Tb,[Point(log1p(ΨmacMat[iiΨmacMat]), 0.0), Point(log1p(ΨmacMat[iiΨmacMat]), 1.0)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
	
				end # for(iiΨmacMat, iΨmacMat) in enumerate(ΨmacMat)

				Legend(Fig[iiσ,N_ΨmacMat+1], Axis_KunsatMat_Tb, framecolor=(:grey, 0.5), labelsize=labelsize, valign=:top, padding=5, tellheight=true, tellwidt=true, nbanks=1, backgroundcolor=:seashell
)
			end # for (iiσ, iσ) in enumerate(σ)

			# General
				resize_to_layout!(Fig)
				trim!(Fig.layout)
				colgap!(Fig.layout, 10)
				rowgap!(Fig.layout, 10)

				Path = raw"D:\TEMP\Plots\MacroKunsat.svg"
				save(Path, Fig)
				display(Fig)
		
	return nothing
	end  # function: name
	# ------------------------------------------------------------------


end #module pumiceManuscript
# ------------------------------------------------------------------


 pumiceManuscript.PLOTTING_KUNSAT_MACRO()

#   include(raw"D:\MAIN\MODELS\AquaPore_Toolkit\src\Temporary\Manuscript\PumiceManuscript.jl")



# include( raw"D:\MAIN\MODELS\AquaPore_Toolkit\src\Temporary\Manuscript\PumiceManuscript.jl")

		
		# Axis_∂θ∂R = Axis(Fig[1,1], xlabel="R [mm]", ylabel="∂θ∂R", xscale=log,  xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5))

		# # , yscale=log10, titlesize=30, xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=true, xticklabelrotation = pi/2, 
		# 	# Axis_∂θ∂R.xticks = (cst.Y ./ TableComplete_θΨ, string.((round.(cst.Y ./ TableComplete_θΨ, digits = 2))))

		# 	# xlims!(Axis_∂θ∂R, cst.Y / Ψ_Max, cst.Y / Ψ_Min)
		# 	Plot1=lines!(Axis_∂θ∂R, cst.Y ./  Ψ , ∂θ∂R_1,  color=:green, linewidth=Linewidth, label="")
		# 	Plot2=lines!(Axis_∂θ∂R,  cst.Y ./ Ψ , ∂θ∂R_2, color=:blue, linewidth=Linewidth, label="")

		# Axis_∂θ∂Ψ = Axis(Fig[1,1], xlabel="Ψ [mm]", ylabel="∂θ∂Ψ", xscale=log)
		# 	# titlesize=30, xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false,xminorticksvisible=true, yminorticksvisible=true, xticklabelrotation = pi/2
			
		# 	# Axis_∂θ∂Ψ.xticks = (TableComplete_θΨ, string.(round.(TableComplete_θΨ, digits=1)))
		# 	# xlims!(Axis_∂θ∂Ψ, Ψ_Min, Ψ_Max)
		# 	Plot3 = lines!(Axis_∂θ∂Ψ, Ψ , ∂θ∂Ψ_1,  color=:green, linewidth=Linewidth, label="")
		# 	Plot4 = lines!(Axis_∂θ∂Ψ, Ψ , ∂θ∂Ψ_2, color=:blue, linewidth=Linewidth, label="")

			# lines!(Axis_∂θ∂Ψ, [Point(Ψm_Min_Mode,0), Point(Ψm_Min_Mode, 1)], color=:green, linewidth=Linewidth)
			# lines!(Axis_∂θ∂Ψ, [Point(Ψm_Max_Mode,0), Point(Ψm_Max_Mode, 1)], color=:blue, linewidth=Linewidth)
			# lines!(Axis_∂θ∂Ψ, [Point(ΨmMac_Mode,0), Point(ΨmMac_Mode, 1)], color=:navyblue, linewidth=Linewidth)
			# lines!(Axis_∂θ∂Ψ, [Point(ΨmacMat,0), Point(ΨmacMat, 0.5)], color=:red, linewidth=Linewidth)

			# lines!(Axis_∂θ∂Ψ, [Point(ΨmMac_Pσ₂_Plus,0), Point(ΨmMac_Pσ₂_Plus, 1.0)], color=:violet)
			# lines!(Axis_∂θ∂Ψ, [Point(ΨmMac_Pσ₂_Minus,0), Point(ΨmMac_Pσ₂_Minus, 1.0)], color=:violet)

			# lines!(Axis_∂θ∂Ψ, [Point(ΨmMac_Pσ₃_Plus,0), Point(ΨmMac_Pσ₃_Plus, 0.8)], color=:violet)
			# lines!(Axis_∂θ∂Ψ, [Point(ΨmMac_Pσ₃_Minus,0), Point(ΨmMac_Pσ₃_Minus, 0.8)], color=:violet)


			
			# lines!(Axis_∂θ∂Ψ, [Point(Ψm_Pσ₂_Plus,0), Point(Ψm_Pσ₂_Plus, 1.0)], color=:blue)
			# lines!(Axis_∂θ∂Ψ, [Point(Ψm_Pσ₂_Minus,0), Point(Ψm_Pσ₂_Minus, 1.0)], color=:blue)

			# lines!(Axis_∂θ∂Ψ, [Point(Ψm_Pσ₃_Plus,0), Point(Ψm_Pσ₃_Plus, 0.8)], color=:blue)
			# lines!(Axis_∂θ∂Ψ, [Point(Ψm_Pσ₃_Minus,0), Point(Ψm_Pσ₃_Minus, 0.8)], color=:blue)

						# Modes
			# Ψm_Min_Mode = hydroRelation.FUNC_ΨmMode(;Ψm₀=Ψm_Min, σ₀=σ)
			# Ψm_Max_Mode = hydroRelation.FUNC_ΨmMode(;Ψm₀=Ψm_Max, σ₀=σ)
			# ΨmMac_Mode  = hydroRelation.FUNC_ΨmMode(;Ψm₀=ΨmMac, σ₀=σMac)
			# Ψm_Mode     = hydroRelation.FUNC_ΨmMode(;Ψm₀=Ψm_Max, σ₀=σ)

		# For plotting
			# ΨmMac_Pσ₂_Plus  = exp(log(ΨmMac_Mode) + Pσ₂ * σMac)
			# ΨmMac_Pσ₂_Minus = exp(log(ΨmMac_Mode) - Pσ₂ * σMac)

			# ΨmMac_Pσ₃_Plus  = exp(log(ΨmMac_Mode) + Pσ₃ * σMac)
			# ΨmMac_Pσ₃_Minus = exp(log(ΨmMac_Mode) - Pσ₃ * σMac)

			# Ψm_Pσ₂_Plus     = exp(log(Ψm_Mode) + Pσ₂ * σ)
			# Ψm_Pσ₂_Minus    = exp(log(Ψm_Mode) - Pσ₂ * σ)

			# Ψm_Pσ₃_Plus     = exp(log(Ψm_Mode) + Pσ₃ * σ)
			# Ψm_Pσ₃_Minus    = exp(log(Ψm_Mode) - Pσ₃ * σ)




# 				using CairoMakie, LaTeXStrings, SpecialFunctions

# x = 0.1:0.1:15

