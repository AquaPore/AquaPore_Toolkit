	using CairoMakie, ColorSchemes
	import SpecialFunctions: erfc, erfcinv


	Path =pwd()
	Path = Path * "//src//"

	include(Path * "Cst.jl")
	include(Path * "hydro\\Wrc.jl")
	include(Path * "hydro\\Kunsat.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : PLOTTING_PORESIZE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function PLOTTING_PORESIZE()

	Pσ₂ = 2.0
	Pσ₃ = 3.0
	Pσ₄ = 4.0

	θs = 0.4
	θr = 0.1
	σ = 1.0
	ΨmacMat = 100.0
	θsMacMat_η = 0.75 
	Ks = 0.008

	θsMacMat = (θs - θr) * θsMacMat_η + θr

	Ψmₐ = √ΨmacMat * exp(σ * Pσ₄)
	Ψmᵦ = ΨmacMat * exp(σ * Pσ₄)

	Ψmₐ_Mode = exp(log(Ψmₐ) - σ^2)
	Ψmᵦ_Mode = exp(log(Ψmᵦ) - σ^2)

	ΨmMac = exp(log(ΨmacMat) / 2.0)

	σMac = log(ΨmacMat) / (2.0 * Pσ₃)

	Rmₐ = cst.Y / Ψmₐ
	Rmᵦ = cst.Y / Ψmᵦ

	Rmₐ_Mode = cst.Y / Ψmₐ_Mode
	Rmᵦ_Mode = cst.Y / Ψmᵦ_Mode

	RmMac = cst.Y / ΨmMac

# -----------------------

	Ψ = 10.0.^(collect(-1:0.000001:6))
	N = length(Ψ)
	∂θ∂Ψ_1 = fill(0.0::Float64 , N)
	∂θ∂Ψ_2 = fill(0.0::Float64 , N)

	Ψ_2_θDual_1 = fill(0.0::Float64 , N)
	Ψ_2_θDual_2 = fill(0.0::Float64 , N)

	Ψ_2_KUNSAT_1 = fill(0.0::Float64 , N)
	Ψ_2_KUNSAT_2 = fill(0.0::Float64 , N)

	∂θ∂R_1 = fill(0.0::Float64 , N) 
	∂θ∂R_2 = fill(0.0::Float64 , N)  

	for iΨ=1:N

		∂θ∂R_1[iΨ] = -wrc.kg.∂θ∂R(R₁= cst.Y/Ψ[iΨ], θs=θs, θr=θr, Rm=Rmₐ , σ=σ, θsMacMat=θsMacMat, RmMac=RmMac, σMac=σMac)

		∂θ∂R_2[iΨ] = -wrc.kg.∂θ∂R(R₁=cst.Y/Ψ[iΨ], θs=θs, θr=θr,  Rm=Rmᵦ, σ=σ, θsMacMat=θsMacMat, RmMac=RmMac, σMac=σMac)

		∂θ∂Ψ_1[iΨ] = wrc.kg.∂θ∂Ψ_NORM(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψmₐ, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, σMac=σMac)

		∂θ∂Ψ_2[iΨ] =  wrc.kg.∂θ∂Ψ_NORM(Ψ₁=Ψ[iΨ], θs=θs, θr=θr,  Ψm=Ψmᵦ, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, σMac=σMac)

		Ψ_2_θDual_1[iΨ] = wrc.kg.Ψ_2_θ(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψmₐ, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, σMac=σMac)
		
		Ψ_2_θDual_2[iΨ] = wrc.kg.Ψ_2_θ(Ψ₁=Ψ[iΨ], θs=θs, θr=θr,  Ψm=Ψmᵦ, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, σMac=σMac)

		Ψ_2_KUNSAT_1[iΨ] = kunsat.kg.Ψ_2_KUNSAT(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψmₐ, σ=σ, θsMacMat=θsMacMat, ΨmMac=ΨmMac, σMac=σMac, Ks=Ks)

		Ψ_2_KUNSAT_2[iΨ] = kunsat.kg.Ψ_2_KUNSAT(Ψ₁=Ψ[iΨ], θs=θs, θr=θr, Ψm=Ψmᵦ, σ=σ, θsMacMat=θsMacMat,ΨmMac=ΨmMac, σMac=σMac, Ks=Ks)
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

	CairoMakie.activate!(type = "svg")
	Fig = Figure(resolution = (5000, 4000), font="Sans", fontsize=NumberSize)
	
	# TableComplete_θΨ = [0.001, 1.0, 10.0, 50, 100.0, 250, 500.0, 1000.0, 2500, 5000.0,100_00.0, 150_00.0,250_00.0, 500_00.0, 1000_00.0,  1500_00.0, 2000_00.0]
	
	Axis0 = Axis(Fig[1,1], xlabel="R [mm]", ylabel="∂θ∂R", xscale=log,  xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5))

	# , yscale=log10, titlesize=30, xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=true, xticklabelrotation = pi/2, 
		# Axis0.xticks = (cst.Y ./ TableComplete_θΨ, string.((round.(cst.Y ./ TableComplete_θΨ, digits = 2))))

		# xlims!(Axis0, cst.Y / Ψ_Max, cst.Y / Ψ_Min)
		Plot1=lines!(Axis0, cst.Y ./  Ψ , ∂θ∂R_1,  color=:green, linewidth=Linewidth, label="")
		Plot2=lines!(Axis0,  cst.Y ./ Ψ , ∂θ∂R_2, color=:blue, linewidth=Linewidth, label="")

	Axis1 = Axis(Fig[2,1], xlabel="Ψ [mm]", ylabel="∂θ∂Ψ", xscale=log)
		# titlesize=30, xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false,xminorticksvisible=true, yminorticksvisible=true, xticklabelrotation = pi/2
		
		# Axis1.xticks = (TableComplete_θΨ, string.(round.(TableComplete_θΨ, digits=1)))
		# xlims!(Axis1, Ψ_Min, Ψ_Max)
		Plot1=lines!(Axis1, Ψ , ∂θ∂Ψ_1,  color=:green, linewidth=Linewidth, label="")
		Plot2=lines!(Axis1, Ψ , ∂θ∂Ψ_2, color=:blue, linewidth=Linewidth, label="")

	Axis2 = Axis(Fig[3,1], xlabel="Ψ [mm]", ylabel="θ [L³ L⁻³]", xscale=log, xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(5))

	# titlesize=30, xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false, yminorticksvisible=true, xticklabelrotation = pi/2

		# xlims!(Axis2, Ψ_Min, Ψ_Max)
		ylims!(Axis2, 0., 0.5)
		# Axis2.xticks = (TableComplete_θΨ, string.( round.(TableComplete_θΨ, digits = 2)))
		lines!(Axis2, Ψ, Ψ_2_θDual_1, color=:green, linewidth=Linewidth)
		lines!(Axis2,Ψ, Ψ_2_θDual_2, color=:blue, linewidth=Linewidth)

	Axis3 = Axis(Fig[4,1], xlabel="Ψ [mm]", ylabel="K(Ψ) [cm h⁻¹]", xscale=log, xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(5))

	# titlesize=30, xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false, xminorticksvisible=true, yminorticksvisible=true, xticklabelrotation = pi/2
		# Axis3.xticks = (TableComplete_θΨ, string.( round.(TableComplete_θΨ, digits=1)))
		lines!(Axis3, Ψ , cst.MmS_2_CmH .* Ψ_2_KUNSAT_1, color=:green, linewidth=Linewidth)
		lines!(Axis3, Ψ , cst.MmS_2_CmH .* Ψ_2_KUNSAT_2, color=:blue, linewidth=Linewidth)
		# xlims!(Axis3, Ψ_Min, Ψ_Max)

		
		
		Legend(Fig[1, 2], Axis1, valign=:top, padding = (0, 0, 20, 0))
		Ψm1 = Int(floor(Ψmₐ))
		Ψm2 = Int(floor(Ψmᵦ))

		axislegend(Axis1, [Plot1, Plot2], ["σ=$σ; Ψm =$Ψm1 mm", "σ=$σ; Ψm =$Ψm2 Kpa"], position = :rt)

		resize_to_layout!(Fig)
		trim!(Fig.layout)
		colgap!(Fig.layout, 50)
		rowgap!(Fig.layout, 50)

		# Path = raw"D:\TEMP\Plots\DthetaH.svg"
		# save(Path, Fig)

		display(Fig)
	
return
end  # function: PLOTTING_PORESIZE
# ------------------------------------------------------------------


PLOTTING_PORESIZE()