# Path =pwd()
	# Path = Path * "//src//"

Path =  "D:\\MAIN\\MODELS\\AquaPore_Toolkit\\src\\"

include(Path * "Including.jl")
include(Path * "Cst.jl")
include(Path * "hydro//HydroRelation.jl")
include(Path * "hydro//Wrc.jl")
include(Path * "hydro//Kunsat.jl")


module pumiceManuscript
	using CairoMakie, ColorSchemes
	import SpecialFunctions: erfc, erfcinv
	import ...wrc, ...kunsat, ...hydroRelation, ..cst

	export PLOTTING_PORESIZE

	Pσ₁ = 1.0
	Pσ₂ = 2.0
	Pσ₃ = 3.0
	Pσ₄ = 4.0

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : PLOTTING_PORESIZE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function PLOTTING_PORESIZE()

	"Units are in mm and S"

	# Input hydro parameters
		θs         = 0.5
		θr         = 0.0
		σ          = 4.0
		θsMacMat_η = 0.75
		θsMacMat   = (θs - θr) * θsMacMat_η + θr
		Ks         = 0.008
		ΨmacMat    = 100.0

	# Deriving macropore hydraulic parameters from ΨmacMat
		σMac    = hydroRelation.FUNC_ΨmacMat_2_σMac(;ΨmacMat=ΨmacMat, Pσ_Mac=3)
		ΨmMac   = hydroRelation.FUNC_ΨmacMat_2_ΨmMac(;ΨmacMat=ΨmacMat, σMac)

		Ψm_Min  = hydroRelation.FUNC_σ_2_Ψm(;ΨmacMat= √ΨmacMat, σ, Pσ=Pσ₃)
		Ψm_Max  = hydroRelation.FUNC_σ_2_Ψm(;ΨmacMat=ΨmacMat, σ, Pσ=Pσ₃)
			
		# Tortuosity parameters
			τa         = 0.5
			τb         = 0.7
			τc         = 0.65
			τaMac      = 0.5
			τbMac      = 0.87
			τcMac      = 0.
			
			Tb_Max = 3.0; Tc_Max = 4.0
			Tb    = Tb_Max * (1.0 - τb)
			TbMac = Tb_Max * (1.0 - τbMac)
			Tc    = Tc_Max * (1.0 - τc)
			# Tc= σ ^ -0.59
			TcMac = Tc_Max * (1.0 - τcMac)
			
			KsMac_Min, KsMat_Min = kunsat.kg.KS_MATMAC_ΨmacMat(θs, θsMacMat, θr, Ψm_Min, σ, ΨmMac, σMac, Ks, Tb, Tc, TbMac, TcMac,"ΨmacMat")
			KsMac_Max, KsMat_Max = kunsat.kg.KS_MATMAC_ΨmacMat(θs, θsMacMat, θr, Ψm_Max, σ, ΨmMac, σMac, Ks, Tb, Tc, TbMac, TcMac,"ΨmacMat")


			θ_Inlection     = (θsMacMat - θr) * 0.5 + θr

		# filling the data -----------------------

			Ψ_Max_Log = log10(exp(log(Ψm_Max) + 2.0 * σ))
			Ψ_Min_Log = log10(exp(log(ΨmMac) - 3.0 * σMac))

			Ψ_Min_Log = log10(0.0001)
			Ψ_Max_Log =log10(1500_00)

			Ψ = 10.0.^(collect(Ψ_Min_Log:0.0001:Ψ_Max_Log))
			N = length(Ψ)

			∂θ∂Ψ_1       = fill(0.0::Float64 , N)
			∂θ∂Ψ_2       = fill(0.0::Float64 , N)

         θDual_Min  = Array{Float64}(undef, N)
         θDual_Max  = Array{Float64}(undef, N)

         Kunsat_Min = Array{Float64}(undef, N)
         Kunsat_Max = Array{Float64}(undef, N)

			for iΨ=1:N
		
				# ∂θ∂Ψ_1[iΨ] = wrc.kg.∂θ∂Ψ_NORM(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Min, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac)

				# ∂θ∂Ψ_2[iΨ] =  wrc.kg.∂θ∂Ψ_NORM(Ψ₁=Ψ[iΨ], θs=θs, θr=θr,  Ψm=Ψm_Max, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac)

				# latexstring("Y_{$(ν)}(x)")

				θDual_Min[iΨ] = wrc.kg.Ψ_2_θ(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Min, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac, KosugiModel_θΨ⍰="ΨmacMat")
				
				θDual_Max[iΨ] = wrc.kg.Ψ_2_θ(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Max, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac, KosugiModel_θΨ⍰="ΨmacMat")

				Kunsat_Min[iΨ] = kunsat.kg.KUNSAT_θΨSe(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Min, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac, Ks=Ks, τa=τa, τb=τb, τc=τc, τaMac=τaMac, τbMac=τbMac, τcMac=τcMac, Option_KosugiModel_KΨ⍰="ΨmacMat", KosugiModel_θΨ⍰="ΨmacMat")

				Kunsat_Max[iΨ] = kunsat.kg.KUNSAT_θΨSe(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Max, σ=σ, θsMacMat=θsMacMat,ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac, Ks=Ks, τa=τa, τb=τb, τc=τc, τaMac=τaMac, τbMac=τbMac, τcMac=τcMac, Option_KosugiModel_KΨ⍰="ΨmacMat", KosugiModel_θΨ⍰="ΨmacMat")
			end

		# ---------------
				# Dimension of figure
               ColourOption_No    = 1
               Linewidth          = 2
               height             = 200
               labelsize          = 12
               textcolor          = :blue
               textsize           = 17
               titlesize          = 20
               width              = height * 3.0
               xgridvisible       = true
               xlabelSize         = 20
               xminorticksvisible = true
               xtickalign         = 1 # 0 is inside and 1 is outside
               xticklabelrotation = π/4
               xticksize          = 10
               xticksmirrored     = false
               xtickwidt          = 0.5
               xtrimspine         = false
               ygridvisible       = false
               ylabelsize         = xlabelSize
               yminorticksvisible = true
               ytickalign         = xtickalign
               yticksize          = xticksize 
               yticksmirrored     = false
               ytickwidt          = xtickwidt 
               ytrimspine         = false
					xgridstyle = :dash 
					ygridstyle = :dash

					ColourOption = [:glasbey_hv_n256,:okabe_ito,:seaborn_bright,:seaborn_colorblind,:seaborn_dark,:seaborn_deep,:tab10,:tableau_10,:tableau_superfishel_stone,:tol_bright]

					Colormap = cgrad(colorschemes[ColourOption[ColourOption_No]], size(colorschemes[ColourOption[ColourOption_No]]), categorical = true)



					Ψticks = [0, 50, 100, 500, 1000,5000,100_00, 500_00, 1000_00,  2000_00]


		CairoMakie.activate!(type="svg", pt_per_unit = 1)
		Fig =  Figure(figure_padding = 10, title="Macropor model", ; fonts = ( ; regular="CMU Serif")) 
		
		Label(Fig[1, 1, Top()], L"Functions $θ(ψ)$ & $K(ψ)$ macropore", valign =:bottom, font =:bold, padding = (0, 0, 5, 0), fontsize=titlesize)


		Ψ_Log = Array{Float64}(undef, N)
		for iZ=1:N
			Ψ_Log[iZ] = log1p(Ψ[iZ])
		end
		

		Axis_θΨ = Axis(Fig[1,1], xlabel= L"$ψ$ [kPa]", ylabel=L"$θ$ [L³ L⁻³]", xticklabelrotation=xticklabelrotation, ylabelsize=ylabelsize, xlabelsize=xlabelSize, xticksize=xticksize, yticksize=yticksize, width=width, height=height,   titlesize=30,  xgridvisible=xgridvisible, ygridvisible=ygridvisible, xminorticksvisible=xminorticksvisible, yminorticksvisible=yminorticksvisible, xtickwidth=xtickwidt, ytickwidth=ytickwidt, xtickalign=xtickalign, ytickalign=ytickalign, xticksmirrored=xticksmirrored, yticksmirrored=yticksmirrored,  xtrimspine=xtrimspine,  ytrimspine=ytrimspine, xgridstyle=xgridstyle, ygridstyle=ygridstyle, yminorticks=IntervalsBetween(5))

			Axis_θΨ.xticks = (log1p.(Ψticks), string.(cst.Mm_2_kPa .* Ψticks))

			hidexdecorations!(Axis_θΨ, ticks=false, grid=false)

			lines!(Axis_θΨ, Ψ_Log, θDual_Min, linewidth=Linewidth, label="σ =$σ, ψₘ_Min", color=Colormap[1])
			lines!(Axis_θΨ, Ψ_Log, θDual_Max, linewidth=Linewidth, label="σ =$σ, ψₘ_Max",  color=Colormap[2])

			lines!(Axis_θΨ,[Point(log1p(ΨmacMat),0), Point(log1p(ΨmacMat), θs)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
			lines!(Axis_θΨ,[Point(log1p(0), θsMacMat ), Point(log1p(ΨmacMat), θsMacMat)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
			lines!(Axis_θΨ,[Point(log1p(0), θs), Point(log1p(ΨmacMat), θs)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)

			text!(log1p(0), θsMacMat, text =L"θ_{sMacMat}", align=(:left,:top),  color=textcolor, fontsize=textsize)
			text!(log1p(0), θs, text =L"θ_{s}", align=(:left,:top),  color=textcolor, fontsize=textsize)
			text!(log1p(ΨmacMat), 0, text =L"ψ_{macMat}", align=(:left,:bottom), rotation = π/2,  color=textcolor, fontsize=textsize)

			Legend(Fig[3,1], Axis_θΨ, framecolor=(:grey, 0.5), labelsize=labelsize, valign=:top, padding=5, tellheight=true, tellwidt=true, nbanks=2, backgroundcolor = (:grey90, 0.25))

			
			Label(Fig[1, 1, TopRight()], "(A)", fontsize=18, font=:bold, padding=(-50, 5, -50, 5), halign=:right)


		Axis_KΨ = Axis(Fig[2,1], xlabel=L"$ψ$ [kPa]", ylabel=L"$K(ψ)$ [cm h⁻¹]", xticklabelrotation=xticklabelrotation, ylabelsize=ylabelsize,	xlabelsize=xlabelSize, xticksize=xticksize, yticksize=yticksize, width=width, height=height, xgridvisible=xgridvisible, ygridvisible=ygridvisible, xminorticksvisible=xminorticksvisible, yminorticksvisible=yminorticksvisible, xtickwidth=xtickwidt, ytickwidth=ytickwidt, xticksmirrored=xticksmirrored,  yticksmirrored=yticksmirrored, xtrimspine=xtrimspine,  ytrimspine=ytrimspine, xtickalign=xtickalign, ytickalign=ytickalign, yminorticks=IntervalsBetween(5), xgridstyle=xgridstyle, ygridstyle=ygridstyle, )

			Axis_KΨ.xticks = (log1p.(Ψticks), string.( cst.Mm_2_kPa .* Ψticks))

			lines!(Axis_KΨ, Ψ_Log , cst.MmS_2_CmH .* Kunsat_Min,  linewidth=Linewidth, color=Colormap[1])
			lines!(Axis_KΨ, Ψ_Log , cst.MmS_2_CmH .* Kunsat_Max, linewidth=Linewidth, color=Colormap[2])

			lines!(Axis_KΨ,[Point(log1p(ΨmacMat),0), Point(log1p(ΨmacMat), cst.MmS_2_CmH .* Ks)], color=:navyblue, linewidth=Linewidth/2.0,  linestyle=:dash)
			lines!(Axis_KΨ,[Point(log1p(0), cst.MmS_2_CmH .* KsMat_Min ), Point(log1p(ΨmacMat), cst.MmS_2_CmH .*  KsMat_Min)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
			lines!(Axis_KΨ,[Point(log1p(0), cst.MmS_2_CmH .* KsMat_Max ), Point(log1p(ΨmacMat), cst.MmS_2_CmH .*  KsMat_Max)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)
			lines!(Axis_KΨ,[Point(log1p(0),  cst.MmS_2_CmH .* Ks), Point(log1p(ΨmacMat), cst.MmS_2_CmH .* Ks)], color=:navyblue, linewidth=Linewidth/2.0, linestyle=:dash)

			text!(log1p(0),  cst.MmS_2_CmH .* KsMat_Min, text =L"K_{sMacMin}", align=(:left,:top), color=textcolor, fontsize=textsize)
			text!(log1p(0),  cst.MmS_2_CmH .* KsMat_Max, text =L"K_{sMacMax}", align=(:left,:top), color=textcolor, fontsize=textsize)
			text!(log1p(0),  cst.MmS_2_CmH .* Ks, text =L"K_{s}", align=(:left,:top), color=textcolor, fontsize=textsize)
			text!(log1p(ΨmacMat), 0, text =L"ψ_{macMat}", align=(:left,:bottom), rotation = π/2, color=textcolor, fontsize=textsize)

			Label(Fig[2, 1, TopRight()], "(B)", fontsize=18, font=:bold, padding = (-50, 5, -50, 5), halign=:right)



			# General
				resize_to_layout!(Fig)
				trim!(Fig.layout)
				colgap!(Fig.layout, 10)
				rowgap!(Fig.layout, 10)

				Path = raw"D:\TEMP\Plots\MacroporeThetaH.svg"
				save(Path, Fig)
				display(Fig)

		
	return
	end  # function: PLOTTING_PORESIZE
end #module pumiceManuscript
# ------------------------------------------------------------------


pumiceManuscript.PLOTTING_PORESIZE()

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

