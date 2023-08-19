	using CairoMakie, ColorSchemes
	import SpecialFunctions: erfc, erfcinv

	include("D:\\MAIN\\MODELS\\AquaPore_Toolkit\\src\\Cst.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : PLOTTING_PORESIZE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function PLOTTING_PORESIZE()

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂θ∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂θ∂ΨMODEL(Ψ₁, θs, θr, Ψm, σ, θsMacMat, ΨmMac, σMac)

			Ψmod_Mat = exp(log(Ψm) - σ^2)

			Ψmod_Mac = exp(log(ΨmMac) - σMac^2)

			# Ψ₁ = max(eps(10.0), Ψ₁)
		
				# If Ψ₁ is positive than ∂θ∂Ψ_Mat should be positive
				∂θ∂Ψ_Mat = (θsMacMat - θr) * exp( -((log(Ψ₁ / Ψm)) ^ 2.0) / (2.0 * σ^2.0)) / (Ψ₁ * σ * √(π * 2.0))

				∂θ∂Ψ_Mat_Mod = (θsMacMat - θr) * exp( -((log(Ψmod_Mat / Ψm)) ^ 2.0) / (2.0 * σ^2.0)) / (Ψmod_Mat * σ * √(π * 2.0))

				∂θ∂Ψ_Mac = (θs - θsMacMat) * exp( -((log(Ψ₁ / ΨmMac)) ^ 2.0) / (2.0 * σMac^2.0)) / (Ψ₁ * σMac * √(π * 2.0))

				∂θ∂Ψ_Mac_Mod = (θs - θsMacMat) * exp( -((log(Ψmod_Mac / ΨmMac)) ^ 2.0) / (2.0 * σMac^2.0)) / (Ψmod_Mac * σMac * √(π * 2.0))

	
		return (∂θ∂Ψ_Mat / ∂θ∂Ψ_Mat_Mod + ∂θ∂Ψ_Mac / ∂θ∂Ψ_Mac_Mod) * 0.5
		end # function ∂θ∂Ψ
	#-------------------------------------------------------------------

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_θDual
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_θDualMODEL(Ψ₁, θs, θr, Ψm, σ, θsMacMat, ΨmMac, σMac)

			θ_Mat = 0.5 * (θsMacMat - θr) * erfc((log( Ψ₁ / Ψm)) / (σ * √2.0)) + θr

			θ_Mac = 0.5 * (θs - θsMacMat) * erfc((log(Ψ₁ / ΨmMac)) / (σMac * √2.0))
		return θ_Mac + θ_Mat
		end # function Ψ_2_θDual
	#-----------------------------------------------------------------

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_KUNSAT_MODEL(Ψ₁, θs, θr, Ψm, σ, θsMacMat, ΨmMac, σMac, Ks)

			θ = Ψ_2_θDualMODEL(Ψ₁, θs, θr, Ψm, σ, θsMacMat, ΨmMac, σMac)

			Se = (θ - θr) / (θs - θr)

			KsMat = Ks * min(max((θsMacMat - θr) / (θs - θr), 0.0), 1.0)			
			Kunsat_Mat =  KsMat * √Se * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2.0

			KsMac = max(Ks - KsMat, 0.0)
			Kunsat_Mac =  KsMac * √Se * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2.0

		return Kunsat_Mat + Kunsat_Mac
		end # function Ψ_2_KUNSAT
	#-------------------------------------------------------------------

	Pσ = 4.0
	ΨmacMat = 100.0
	ΨmMac = √ΨmacMat
	σMac = log(ΨmacMat) / (2.0 * 3)
	θs = 0.5
	θr = 0
	σ = 1.5
	Ψm = (√ΨmacMat + ΨmacMat) * 0.5 * exp(σ * Pσ)
	θsMacMat = θs * 0.75
	Ks = 0.008

# -----------------------

	Ψ = 10.0.^(collect(0:0.001:6))
	N = length(Ψ)
	∂θ∂Ψ_1 = fill(0.0::Float64 , N)
	∂θ∂Ψ_2 = fill(0.0::Float64 , N)

	Ψ_2_θDual_1= fill(0.0::Float64 , N)
	Ψ_2_θDual_2= fill(0.0::Float64 , N)

	Ψ_2_KUNSAT_1 = fill(0.0::Float64 , N)
	Ψ_2_KUNSAT_2 = fill(0.0::Float64 , N)  

	for iΨ=1:N
		∂θ∂Ψ_1[iΨ] = ∂θ∂ΨMODEL(Ψ[iΨ], θs, θr, ΨmacMat * exp(σ * Pσ), σ, θsMacMat, ΨmMac, σMac)
		∂θ∂Ψ_2[iΨ] = ∂θ∂ΨMODEL(Ψ[iΨ], θs, θr, √ΨmacMat * exp(σ * Pσ), σ, θsMacMat, ΨmMac, σMac)

		Ψ_2_θDual_1[iΨ] = Ψ_2_θDualMODEL(Ψ[iΨ], θs, θr, ΨmacMat * exp(σ * Pσ), σ, θsMacMat, ΨmMac, σMac)
		Ψ_2_θDual_2[iΨ] = Ψ_2_θDualMODEL(Ψ[iΨ], θs, θr, √ΨmacMat * exp(σ * Pσ), σ, θsMacMat, ΨmMac, σMac)

		Ψ_2_KUNSAT_1[iΨ] = Ψ_2_KUNSAT_MODEL(Ψ[iΨ], θs, θr, ΨmacMat * exp(σ * Pσ), σ, θsMacMat, ΨmMac, σMac, Ks)

		Ψ_2_KUNSAT_2[iΨ] = Ψ_2_KUNSAT_MODEL(Ψ[iΨ], θs, θr, √ΨmacMat * exp(σ * Pσ), σ, θsMacMat, ΨmMac, σMac, Ks)
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
        

	CairoMakie.activate!(type = "svg")
	Fig = Figure(resolution = (5000, 4000), font="Sans", fontsize=NumberSize)
	
	TableComplete_θΨ = [0.0, 1.0, 10.0, 50, 100.0, 250, 500.0, 1000.0, 2500, 5000.0,100_00.0, 250_00.0, 500_00.0, 1000_00.0,  2000_00.0]

	Axis1 = Axis(Fig[1,1], xlabel="Ψ [Kpa]", ylabel="Normalised ∂θ∂Ψ", xscale=Makie.log, titlesize=30, xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false,xminorticksvisible=true, yminorticksvisible=true)
		Axis1.xticks = (TableComplete_θΨ, string.( Int.(TableComplete_θΨ)))
		
		Plot1=lines!(Axis1, cst.Mm_2_kPa * Ψ , ∂θ∂Ψ_1,  color=:green, linewidth=Linewidth, label="")
		Plot2=lines!(Axis1, cst.Mm_2_kPa * Ψ , ∂θ∂Ψ_2, color=:blue, linewidth=Linewidth, label="")

	Axis2 = Axis(Fig[2,1], xlabel="Ψ [kPa]", ylabel="θ(Ψ)", xscale=Makie.log, titlesize=30, xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false,xminorticksvisible=true, yminorticksvisible=true)
		Axis2.xticks = (TableComplete_θΨ, string.( Int.(TableComplete_θΨ)))
		lines!(Axis2, cst.Mm_2_kPa .* Ψ , Ψ_2_θDual_1, color=:green, linewidth=Linewidth)
		lines!(Axis2, cst.Mm_2_kPa .* Ψ , Ψ_2_θDual_2, color=:blue, linewidth=Linewidth)

	Axis3 = Axis(Fig[3,1], xlabel="Ψ [kPa]", ylabel="K(Ψ)", xscale=Makie.log, titlesize=30, xlabelsize=XlabelSize, ylabelsize=YlabelSize, xgridvisible=false, ygridvisible=false, xminorticksvisible=true, yminorticksvisible=true)
		Axis3.xticks = (TableComplete_θΨ, string.( Int.(TableComplete_θΨ)))
		lines!(Axis3, cst.Mm_2_kPa .* Ψ , Ψ_2_KUNSAT_1, color=:green, linewidth=Linewidth)
		lines!(Axis3, cst.Mm_2_kPa .* Ψ , Ψ_2_KUNSAT_2, color=:blue, linewidth=Linewidth)

		Legend(Fig[1, 2], Axis1, valign=:top, padding = (0, 0, 20, 0))

		Ψm1 = Int(floor(cst.Mm_2_kPa*ΨmacMat * exp(σ * Pσ)))
		Ψm2 = Int(floor(cst.Mm_2_kPa*√ΨmacMat * exp(σ * Pσ)))

		axislegend(Axis1, [Plot1, Plot2], ["σ=$σ; Ψm =$Ψm1 Kpa", "σ=$σ; Ψm =$Ψm2 Kpa"], position = :rt)

		resize_to_layout!(Fig)
		trim!(Fig.layout)
		colgap!(Fig.layout, 50)
		rowgap!(Fig.layout, 50)

		Path = raw"D:\TEMPORARY\Plots\DthetaH.svg"
		save(Path, Fig)

	display(Fig)
	
return
end  # function: PLOTTING_PORESIZE
# ------------------------------------------------------------------


PLOTTING_PORESIZE()