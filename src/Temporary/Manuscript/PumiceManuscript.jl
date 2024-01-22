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
			import ..wrc, ..kunsat, ...hydroRelation, ..cst

			export PLOTTING_PORESIZE
			Pσ₁ = 1.0
			Pσ₂ = 2.0
			Pσ₃ = 3.0
			Pσ₄ = 4.0

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOTTING_PORESIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOTTING_PORESIZE()
			"""Units are in mm and S"""

			# Input hydro parameters
            θs         = 0.4
            θr         = 0.0
            σ          = 2.0
            θsMacMat_η = 0.75
            θsMacMat   = (θs - θr) * θsMacMat_η + θr
            Ks         = 0.008
            ΨmacMat    = 100.0

			# Deriving macropore hydraulic parameters from ΨmacMat
				σMac    = hydroRelation.FUNC_ΨmacMat_2_σMac(;ΨmacMat=ΨmacMat, Pσ_Mac=3)

				ΨmMac   = hydroRelation.FUNC_ΨmacMat_2_ΨmMac(;ΨmacMat=ΨmacMat, σMac)

				Ψm_Min  = hydroRelation.FUNC_σ_2_Ψm(;ΨmacMat= ΨmacMat, σ, Pσ=Pσ₃)
			
				Ψm_Max  = hydroRelation.FUNC_σ_2_Ψm(;ΨmacMat=ΨmacMat, σ, Pσ=Pσ₃)
				
			# Tortuosity parameters
            τa         = 0.5
            τb         = 0.5
            τc         = 0.5
            τaMac      = 0.5
            τbMac      = 0.5
            τcMac      = 0.5
				
				Tb_Max = 3.0; Tc_Max = 4.0
				Tb    = Tb_Max * (1.0 - τb)
				TbMac = Tb_Max * (1.0 - τbMac)
				Tc    = Tc_Max * (1.0 - τc)
				# Tc= σ ^ -0.59
				TcMac = Tc_Max * (1.0 - τcMac)

				


				Ks_Mac, Ks_Mat = kunsat.kg.KS_MATMAC_ΨmacMat(θs::Float64, θsMacMat::Float64, θr::Float64, Ψm_Min::Float64, σ::Float64, ΨmMac::Float64, σMac::Float64, Ks::Float64, Tb::Float64, Tc::Float64, TbMac::Float64, TcMac::Float64,"ΨmacMat")


			

            θ_Inlection     = (θsMacMat - θr) * 0.5 + θr

			# filling the data -----------------------

				Ψ_Max_Log = log10(exp(log(Ψm_Max) + 2.0 * σ))
				Ψ_Min_Log = log10(exp(log(ΨmMac) - 3.0 * σMac))
				

				Ψ = 10.0.^(collect(Ψ_Min_Log:0.0001:Ψ_Max_Log))
				N = length(Ψ)

            ∂θ∂Ψ_1       = fill(0.0::Float64 , N)
            ∂θ∂Ψ_2       = fill(0.0::Float64 , N)

            Ψ_2_θDual_1  = fill(0.0::Float64 , N)
            Ψ_2_θDual_2  = fill(0.0::Float64 , N)

            Ψ_2_KUNSAT_1 = fill(0.0::Float64 , N)
            Ψ_2_KUNSAT_2 = fill(0.0::Float64 , N)

				for iΨ=1:N
			
					# ∂θ∂Ψ_1[iΨ] = wrc.kg.∂θ∂Ψ_NORM(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Min, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac)

					# ∂θ∂Ψ_2[iΨ] =  wrc.kg.∂θ∂Ψ_NORM(Ψ₁=Ψ[iΨ], θs=θs, θr=θr,  Ψm=Ψm_Max, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac)

					Ψ_2_θDual_1[iΨ] = wrc.kg.Ψ_2_θ(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Min, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac, KosugiModel_θΨ⍰="ΨmacMat")
					
					Ψ_2_θDual_2[iΨ] = wrc.kg.Ψ_2_θ(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Max, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac, KosugiModel_θΨ⍰="ΨmacMat")

					Ψ_2_KUNSAT_1[iΨ] = kunsat.kg.KUNSAT_θΨSe(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Min, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac, Ks=Ks, τa=τa, τb=τb, τc=τc, τaMac=τaMac, τbMac=τbMac, τcMac=τcMac, Option_KosugiModel_KΨ⍰="ΨmacMat", KosugiModel_θΨ⍰="ΨmacMat")

					Ψ_2_KUNSAT_2[iΨ] = kunsat.kg.KUNSAT_θΨSe(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψm_Max, σ=σ, θsMacMat=θsMacMat,ΨmMac=ΨmMac, ΨmacMat=ΨmacMat, σMac=σMac, Ks=Ks, τa=τa, τb=τb, τc=τc, τaMac=τaMac, τbMac=τbMac, τcMac=τcMac, Option_KosugiModel_KΨ⍰="ΨmacMat", KosugiModel_θΨ⍰="ΨmacMat")
				end

			# ---------------
					# Dimension of figure
						Height = 3000 # Height of plot
						Width  = 5500  # Width of plot

					# Size of X and Y label
						XlabelSize = 125
						YlabelSize = 125
						NumberSize = 80

					# Title
						TitleSize = 60

					# Labels size of colourbar
						TickLabelSize = 35
						TickSize      = 20
						LabelSize     = 35

						Linewidth = 12

					# Limit
						Ψ_Min = 0.001
						Ψ_Max = 2000_00.0

			# CairoMakie.activate!(type = "svg")
			# Fig = Figure(size = (5000, 4000), font="Sans", fontsize=NumberSize)
			Fig = Figure( size=(5000, 4000), font="Sans", fontsize=NumberSize)
			



			Axis_θΨ = Axis(Fig[2,1], xlabel="Ψ [mm]", ylabel="θ [L³ L⁻³]", xscale=log, xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(5))

			# titlesize=30, xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false, yminorticksvisible=true, xticklabelrotation = pi/2

				# xlims!(Axis_θΨ, Ψ_Min, Ψ_Max)
				ylims!(Axis_θΨ, 0., 0.5)
				# Axis_θΨ.xticks = (TableComplete_θΨ, string.( round.(TableComplete_θΨ, digits = 2)))
				lines!(Axis_θΨ, Ψ, Ψ_2_θDual_1, color=:green, linewidth=Linewidth)
				lines!(Axis_θΨ,Ψ, Ψ_2_θDual_2, color=:blue, linewidth=Linewidth)
				lines!(Axis_θΨ, Ψ, (Ψ ./ Ψ ) .* θ_Inlection , color=:yellow, linewidth=Linewidth/2.0)
				lines!(Axis_θΨ,[Point(ΨmacMat,0), Point(ΨmacMat,θsMacMat)], color=:brown, linewidth=Linewidth/2.0)
				lines!(Axis_θΨ, Ψ, (Ψ ./ Ψ ) .* θsMacMat , color=:grey, linewidth=Linewidth/2.0)
				lines!(Axis_θΨ, Ψ, (Ψ ./ Ψ ) .* θs , color=:grey, linewidth=Linewidth/2.0)

			Axis_KΨ = Axis(Fig[4,1], xlabel="Ψ [mm]", ylabel="K(Ψ) [cm h⁻¹]", xscale=log, xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(5))

			# titlesize=30, xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false, xminorticksvisible=true, yminorticksvisible=true, xticklabelrotation = pi/2
				# Axis_KΨ.xticks = (TableComplete_θΨ, string.( round.(TableComplete_θΨ, digits=1)))
				lines!(Axis_KΨ, Ψ , cst.MmS_2_CmH .* Ψ_2_KUNSAT_1, color=:green, linewidth=Linewidth)
				lines!(Axis_KΨ, Ψ , cst.MmS_2_CmH .* Ψ_2_KUNSAT_2, color=:blue, linewidth=Linewidth)
				lines!(Axis_KΨ, Ψ, (Ψ ./ Ψ ) .* cst.MmS_2_CmH .* Ks_Mat, color=:grey, linewidth=Linewidth/2.0)
				lines!(Axis_KΨ, Ψ, (Ψ ./ Ψ ) .* cst.MmS_2_CmH .* Ks, color=:grey, linewidth=Linewidth/2.0)
				lines!(Axis_KΨ,[Point(ΨmacMat,0), Point(ΨmacMat,cst.MmS_2_CmH * Ks_Mat)], color=:brown, linewidth=Linewidth/2.0)
				# xlims!(Axis_KΨ, Ψ_Min, Ψ_Max)

				
				
				# Legend(Fig[1, 3], Axis_∂θ∂Ψ, valign=:top, padding = (0, 0, 20, 0))
				Ψm1 = Int(floor(Ψm_Min))
				Ψm2 = Int(floor(Ψm_Max))

				# axislegend(Axis_∂θ∂Ψ, [Axis_∂θ∂Ψ], ["σ=$σ; Ψm =$Ψm1 mm", "σ=$σ; Ψm =$Ψm2 Kpa"], position = :rt)

				resize_to_layout!(Fig)
				trim!(Fig.layout)
				colgap!(Fig.layout, 50)
				rowgap!(Fig.layout, 50)

				# Path = raw"D:\TEMP\Plots\DthetaH.svg"
				display(Fig)
				# save(Path, Fig)

			
		return
		end  # function: PLOTTING_PORESIZE
end #module pumiceManuscript
# ------------------------------------------------------------------


pumiceManuscript.PLOTTING_PORESIZE()

# include( raw"D:\MAIN\MODELS\AquaPore_Toolkit\src\Temporary\Manuscript\PumiceManuscript.jl")

			# TableComplete_θΨ = [0.001, 1.0, 10.0, 50, 100.0, 250, 500.0, 1000.0, 2500, 5000.0,100_00.0, 150_00.0,250_00.0, 500_00.0, 1000_00.0,  1500_00.0, 2000_00.0]
			
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
				# lines!(Axis_∂θ∂Ψ, [Point(ΨmMac_Mode,0), Point(ΨmMac_Mode, 1)], color=:brown, linewidth=Linewidth)
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