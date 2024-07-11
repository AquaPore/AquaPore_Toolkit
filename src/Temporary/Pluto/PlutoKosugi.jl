### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ fcfdad85-3e5d-4642-93a7-7843bd2cf50f
begin
    using Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate(Base.current_project())

	Path =  "D:\\MAIN\\MODELS\\AquaPore_Toolkit\\src\\"
	include(Path * "Cst.jl")
	include(Path * "hydro//HydroRelation.jl")
	include(Path * "hydro//Wrc.jl")
	include(Path * "hydro//Kunsat.jl")

	using CairoMakie, ColorSchemes, GLMakie
	 using WGLMakie
	using CSV, Tables, DataFrames, LaTeXStrings
	import SpecialFunctions: erfc, erfcinv
	import ..wrc, ..kunsat, ..hydroRelation, ..cst;
end

# ╔═╡ 777b58c7-180f-41a9-bfa5-edbbc246c939
html"""<style>
main {
    max-width: 96%;
    margin-left: 1%;
    margin-right: 2% !important;
}
"""

# ╔═╡ 111ff21d-bc50-4c05-933d-db813261a196
begin

	# Observed data
θΨ_LogΨobs = log1p.([0.0, 40.0, 70.0, 100.0, 500.0, 1000.0, 2000.0, 4000.0, 10000.0, 150000.0]) 
θΨ_θobs_1 = [0.672, 0.612, 0.61, 0.609, 0.588, 0.572, 0.559, 0.554, 0.549, 0.163]
θΨ_θobs_2 = [0.733, 0.507, 0.48, 0.463, 0.356, 0.348, 0.34, 0.331, 0.318, 0.036]
θΨ_θobs_3 = [0.658, 0.573, 0.567, 0.561, 0.516, 0.501, 0.488, 0.488, 0.487, 0.081]

KΨ_Kobs_1 = [0.00320, 0.00360, 0.00300, 0.00220] 
KΨ_Kobs_2 = [0.56300, 0.05200, 0.01400, 0.00340]
KΨ_Kobs_3 = [0.00470, 0.00440, 0.00480, 0.00360];

θΨ_Ψobs_Sand = log1p.([0 
70
170
360
420
490
590
660
810
980
1180
1380
1560
1780
2180
2540
2800
])
θΨ_θobs_Sand = [0.492307692
0.4923
0.492
0.484
0.483
0.466
0.391
0.357
0.296
0.254
0.22
0.198
0.183
0.172
0.152
0.143
0.139
]

θΨ_Ψobs_Silt = log1p.([0
10
50
100
200
400
800
1600
3450
6900
20000
50000
100000
150000
])

θΨ_θobs_Silt = [0.452830189
0.409
0.4
0.386
0.36
0.321
0.28
0.258
0.238
0.221
0.192
0.17
0.146
0.127
]

θΨ_Ψobs_Clay = log1p.([0
30
70
150
210
320
430
520
780
1010
1300
1580
10000
160000
4200000
])
θΨ_θobs_Clay = [0.437735849
0.43
0.415
0.398
0.389
0.377
0.37
0.363
0.355
0.347
0.341
0.338
0.308
0.217
0.042
]



KΨ_LogΨobs = log1p.([10.0, 40.0, 70.0, 100.0])  

KΨ_Ψobs_Sand = log1p.([0
20
80
250
360
490
550
620
740
910
])

KΨ_Kobs_Sand = [0.006805556
0.006793981
0.006805556
0.006759259
0.006736111
0.003541667
0.001724537
0.000892361
0.000393519
0.000148148
]

KΨ_Ψobs_Silt = log1p.([0
50
100
200
400
800
1600
])
KΨ_Kobs_Silt = [0.001290509
0.000379977
0.000275
3.69E-05
1.26E-05
2.24E-06
1.29E-06
]

KΨ_Ψobs_Clay = log1p.( [320
400
500
630
790
1000
1260
1580
2000
2510
3160
3980
5010
6310
7940
10000
])
KΨ_Kobs_Clay = [3.04E-06
2.21E-06
1.57E-06
1.10E-06
7.60E-07
5.17E-07
3.48E-07
2.31E-07
1.53E-07
9.95E-08
6.48E-08
4.28E-08
2.78E-08
1.74E-08
1.16E-08
6.94E-09
]
end

# ╔═╡ 119793bf-274a-43b4-9b8c-0aa656adc71d
begin
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θΨMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θΨMODEL(Ψ; KosugiModel_θΨ⍰="ΨmacMat", Pσ_Mac=2.0, θr, θs, θsMacMat, σ, Ψm=1000.0, ΨmacMat_2_σMac_ΨmMac=true, ΨmacMat=100)
	
			θsim = zeros(length(Ψ))
	
			for (iiΨ, iΨ) in enumerate(Ψ)
				θsim[iiΨ] = wrc.kg.Ψ_2_θ(;Ψ₁=iΨ, θs, θsMacMat, θr, Ψm, σ, ΨmacMat, KosugiModel_θΨ⍰, ΨmacMat_2_σMac_ΨmMac, Pσ_Mac)
			end
			
		return θsim
		end  # function: θΨMODEL
	# ------------------------------------------------------------------
	
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KΨMODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KΨMODEL(Ψ; KosugiModel_KΨ⍰="ΨmacMat", KosugiModel_θΨ⍰="ΨmacMat", KosugiModel_σ_2_Tb=false, Ks, Pσ_Mac=2.0, θr, θs, θsMacMat, σ_Max=4.0, σ_Min=0.7, σ, τa=0.5, τaMac=0.5, τb=1.103, τbMac=0.619, τc=1.0, τcMac=2.0, τₚ=3.0, Ψm, ΨmacMat_2_σMac_ΨmMac=true, ΨmacMat)
	
			Ksim = zeros(length(Ψ))
	
			for (iiΨ, iΨ) in enumerate(Ψ)
				Ksim[iiΨ] = kunsat.kg.KUNSAT_θΨSe(;Ψ₁=iΨ, θs, θsMacMat, θr, Ψm, σ, ΨmacMat, Ks, τa, τb, τc, τₚ, τaMac, τbMac, τcMac, σ_Min, σ_Max, KosugiModel_KΨ⍰, KosugiModel_θΨ⍰, KosugiModel_σ_2_Tb, Pσ_Mac, ΨmacMat_2_σMac_ΨmMac)
			end
	
			KsMat = kunsat.kg.KUNSAT_θΨSe(;Ψ₁=ΨmacMat, θs, θsMacMat, θr, Ψm, σ, ΨmacMat, Ks, τa, τb, τc, τₚ, τaMac, τbMac, τcMac, σ_Min, σ_Max, KosugiModel_KΨ⍰, KosugiModel_θΨ⍰, KosugiModel_σ_2_Tb, Pσ_Mac, ΨmacMat_2_σMac_ΨmMac)
			
		return Ksim, KsMat
		end  # function: KΨMODEL
	# ------------------------------------------------------------------
	
	
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : DISTRIBUTION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function DISTRIBUTION(Ψ; KosugiModel_θΨ⍰="ΨmacMat", Pσ_Mac=2.0, θr, θs, θsMacMat, σ, Ψm=1000.0, ΨmacMat_2_σMac_ΨmMac=true, ΨmacMat=100)
	
		ΨmMac = hydroRelation.FUNC_ΨmacMat_2_ΨmMac(;ΨmacMat)
		σMac  = hydroRelation.FUNC_ΨmacMat_2_σMac(;ΨmacMat)
	
		∂θ∂Ψ_Mat_Vect = zeros(length(Ψ))
		∂θ∂Ψ_Mac_Vect = zeros(length(Ψ))
	
		for (iiΨ, iΨ) in enumerate(Ψ)
	
			Ψmod_Mat = exp(log(Ψm) - σ^2)
	
			Ψmod_Mac = exp(log(ΨmMac) - σMac^2)
	
			if iΨ > eps(100.0)
				∂θ∂Ψ_Mat(Ψ₁) = (θsMacMat - θr) * exp( -((log(Ψ₁ / Ψm)) ^ 2.0) / (2.0 * σ^2.0)) / (Ψ₁ * σ * √(π * 2.0))
	
				∂θ∂Ψ_Mat_Mod = (θsMacMat - θr) * exp( -((log(Ψmod_Mat / Ψm)) ^ 2.0) / (2.0 * σ^2.0)) / (Ψmod_Mat * σ * √(π * 2.0))
	
				∂θ∂Ψ_Mac(Ψ₁) = (θs - θsMacMat) * exp( -((log(Ψ₁ / ΨmMac)) ^ 2.0) / (2.0 * σMac^2.0)) / (Ψ₁ * σMac * √(π * 2.0))
	
				∂θ∂Ψ_Mac_Mod = (θs - θsMacMat) * exp( -((log(Ψmod_Mac / ΨmMac)) ^ 2.0) / (2.0 * σMac^2.0)) / (Ψmod_Mac * σMac * √(π * 2.0))
	
				∂θ∂Ψ_Mat_Vect[iiΨ] = ∂θ∂Ψ_Mat(iΨ) / ∂θ∂Ψ_Mat(Ψmod_Mat)
	
				∂θ∂Ψ_Mac_Vect[iiΨ] = ∂θ∂Ψ_Mac(iΨ) / ( ∂θ∂Ψ_Mac(Ψmod_Mac) + eps())
	
			else
				∂θ∂Ψ_Mat_Vect[iiΨ] = 0.0 
				∂θ∂Ψ_Mac_Vect[iiΨ] = 0.0
	
			end # function ∂θ∂Ψ_NORM
	
		end
		
	return ∂θ∂Ψ_Mat_Vect, ∂θ∂Ψ_Mac_Vect
	end  # function: DISTRIBUTION
	# ------------------------------------------------------------------
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : APPEND_HYDRO
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function APPEND_HYDRO(HydroVect)
			Y = []
			for iHydroVect ∈ HydroVect
				append!(Y, iHydroVect)
			end
		return Y
		end  # function: APPEND_HYDRO
	# ------------------------------------------------------------------
	
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FUNC_σ_2_Ψm
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNC_σ_2_Ψm(ΨmacMat, σ; Pσ=3.0)
			# return ΨmModel =(ΨmacMat ^ 0.75) * exp(σ * Pσ)
	
			return exp(log(√ΨmacMat * exp(σ * Pσ)) + log(ΨmacMat * exp(σ * Pσ)))
		end  # function: FUNC_σ_2_Ψm
	# ------------------------------------------------------------------
end

# ╔═╡ c98dca7b-0612-4c2a-804f-b545acac4b34
begin
	#  For every ψ
	Ψ_Min_Log = log10(0.0001); Ψ_Max_Log = log10(1500_00.0)
	Ψ = 10.0.^(collect(Ψ_Min_Log:0.0001:Ψ_Max_Log))
	N_Ψ = length(Ψ)
	
	Ψ_Log = Array{Float64}(undef, N_Ψ)
	for iZ=1:N_Ψ
		Ψ_Log[iZ] = log1p(Ψ[iZ])
	end
end

# ╔═╡ 0e06f203-8f2a-4844-8678-9b2bf9a3f495
begin
	using  Markdown
	# WGLMakie.aPage() # for Franklin, you still need to configure
WGLMakie.activate!()
Makie.inline!(true)  # Make sure to inline plots into Documenter output!
	
	# WGLMakie.activate!(inline=true) 
	 # GLMakie.activate!(inline=false) 
	
	Linewidth  = 4
	xlabelSize = 30
	xticksize  = 10
	xgridvisible = false
	Width = 800
	Height = 200
	
	Fig = Figure(size = (3000, Height*4.0*1.5))


	
	sg = SliderGrid(Fig[4,1],
	    (label="θs", range=0.25:0.005:0.6, startvalue       = 0.5),
	    (label="θsMacMat", range=0.75:0.005:1.0, startvalue = 0.8),
	    (label="θr", range=0.0:0.005:0.25, startvalue       = 0.1),
	    (label="σ", range=0.75:0.005:3.75, startvalue       = 2.0),
	    (label="Log10_Ψm", range=1.0:0.005:10.0, startvalue  = 5.0),
	    (label="Log10_ΨmacMat", range=0.01:0.005:2.0, startvalue = 0.5),
	    (label="Kₛ", range=1.0:0.5:50.0, startvalue      = 10.0),
	    width = 1000, tellheight = true)
	
	    θs       = sg.sliders[1].value
	    θsMacMat = sg.sliders[2].value
	    θr       = sg.sliders[3].value
	    σ        = sg.sliders[4].value
	    Ψm       = sg.sliders[5].value
	    ΨmacMat  = sg.sliders[6].value
	    Ks       = sg.sliders[7].value
	    
	    obs_func = on(σ) do val
	        val > 0 && Makie.set_close_to!(sg.sliders[5], log10(FUNC_σ_2_Ψm(ΨmacMat[], σ[])))
	    end
	    # obs_func = on(θs) do val
	    #     val < θsMacMat[] && Makie.set_close_to!(sg.sliders[1], θsMacMat[])
	    # end
	
	    # Menu
	         Menu1 = Menu(Fig[1,1:2], options = ["ΨmacMat", "Traditional"], default = "ΨmacMat")
	         Fig[1, 2] = vgrid!(
	         Label(Fig, "Model", width = nothing),
	         Menu1,
	         tellheight = false, width = 150)
	
	         Model = Observable{Any}([""])
	         # Model ="Traditional"
	         on(Menu1.selection) do S
	            Model[]= S 
	         end
	         notify(Menu1.selection)
	      
	      # LIFT DISTRIBUTION
	         ∂θ∂Ψ_Mat = lift((θs, θsMacMat, θr, σ, Ψm, ΨmacMat, Model) -> DISTRIBUTION(Ψ; θs=θs, θsMacMat=θsMacMat*θs, θr=θr, σ=σ, Ψm=10.0^Ψm, ΨmacMat=10.0^ΨmacMat, KosugiModel_θΨ⍰=Model)[1], θs, θsMacMat, θr, σ, Ψm, ΨmacMat, Model)
	
	         ∂θ∂Ψ_Mac = lift((θs, θsMacMat, θr, σ, Ψm, ΨmacMat, Model) -> DISTRIBUTION(Ψ; θs=θs, θsMacMat=θsMacMat*θs, θr=θr, σ=σ, Ψm=10.0^Ψm, ΨmacMat=10.0^ΨmacMat, KosugiModel_θΨ⍰=Model)[2], θs, θsMacMat, θr, σ, Ψm, ΨmacMat, Model)
	   
	      # LIFT θ(ψ)
	        θsim = lift((θs, θsMacMat, θr, σ, Ψm, ΨmacMat, Model) -> θΨMODEL(Ψ; θs=θs, θsMacMat=θsMacMat*θs, θr=θr, σ=σ, Ψm=10.0^Ψm, ΨmacMat=10.0^ΨmacMat, KosugiModel_θΨ⍰=Model), θs, θsMacMat, θr, σ, Ψm, ΨmacMat, Model)
	
	        Y_Line_θΨ = lift((θs, θsMacMat, θr) -> APPEND_HYDRO([θs, θsMacMat*θs, θr]), θs, θsMacMat, θr)
	
	        X_Line_θΨ = lift((ΨmacMat) -> APPEND_HYDRO([log1p(10.0^ΨmacMat)]), ΨmacMat)
	
	    # LIFT K(ψ)
	        Ksim = lift((θs, θsMacMat, θr, σ, Ψm, ΨmacMat, Ks, Model) -> KΨMODEL(Ψ; θs=θs, θsMacMat=θsMacMat*θs, θr=θr, σ=σ, Ψm=10.0^Ψm, ΨmacMat=10.0^ΨmacMat, Ks=Ks, KosugiModel_KΨ⍰=Model, KosugiModel_θΨ⍰=Model)[1], θs, θsMacMat, θr, σ, Ψm, ΨmacMat, Ks, Model)
	
	        KsMat = lift((θs, θsMacMat, θr, σ, Ψm, ΨmacMat, Ks, Model) -> KΨMODEL(Ψ; θs=θs, θsMacMat=θsMacMat*θs, θr=θr, σ=σ, Ψm=10.0^Ψm, ΨmacMat=10.0^ΨmacMat, Ks=Ks, KosugiModel_KΨ⍰=Model, KosugiModel_θΨ⍰=Model)[2], θs, θsMacMat, θr, σ, Ψm, ΨmacMat, Ks, Model)
	
	        Y_Line_KsMat = lift((Ks, KsMat) -> APPEND_HYDRO([Ks, KsMat]), Ks, KsMat)
	
	
	    # PLOTTING AX_1: θ(Ψ)
	        Ax_1 = Axis(Fig[2, 1], width=Width, height=Height, xticklabelrotation = π/4.0, xlabel= L"$ψ$ [kPa]", ylabel=L"$\theta(\psi)$ [L³ L⁻³]",xlabelsize=xlabelSize, ylabelsize=xlabelSize, xticksize=xticksize, xgridvisible=xgridvisible, ygridvisible=xgridvisible)
	            ylims!(Ax_1, 0.0, 0.6)
	        
	            Ψticks = [0, 50, 100, 500, 1000,5000,100_00, 500_00, 1000_00, 1500_00] # mm
	            Ax_1.xticks = (log1p.(Ψticks), string.(cst.Mm_2_kPa .* Ψticks))
	            
	            θticks = [0.0, 0.25, 0.5, 0.75, 1.0] # mm
	            Ax_1.yticks = (θticks, string.(θticks))
	        
	            lines!(Ax_1, Ψ_Log, θsim, linewidth=Linewidth, color=:red2, label="θΨ_MacMat")

	
	            hlines!(Ax_1, Y_Line_θΨ; xmin=0.0, color=:blue2, linewidth=Linewidth/2.0, linestyle=:dash)
	
	            vlines!(Ax_1, X_Line_θΨ; ymin = 0.0, color=:blue2, linewidth=Linewidth/2.0, linestyle=:dash)
	
	            # scatter!(Ax_1, θΨ_LogΨobs, θΨ_θobs_1, marker=:hexagon, markersize=25)
	
	            # scatter!(Ax_1, θΨ_LogΨobs, θΨ_θobs_2, marker=:diamond, markersize=25)
	
	            scatter!(Ax_1, θΨ_Ψobs_Sand, θΨ_θobs_Sand, marker=:diamond, markersize=30, color=:darkgoldenrod3)
	
	            scatter!(Ax_1, θΨ_Ψobs_Silt, θΨ_θobs_Silt, markersize=30, color=:teal)
	
	            # scatter!(Ax_1, θΨ_Ψobs_Clay, θΨ_θobs_Clay, marker=:diamond, markersize=25)
	
	    # PLOTTING AX_2: K(Ψ)
	        Ax_2 = Axis(Fig[3, 1], width=Width, height=Height, xticklabelrotation = π / 4.0, xlabel= L"$ψ$ [kPa]", ylabel=L"$K(\psi)$ [mm hr ⁻¹]", xlabelsize=xlabelSize, ylabelsize=xlabelSize, xgridvisible=xgridvisible, ygridvisible=xgridvisible, yscale=log10)
	
	            Ψticks = [0, 50, 100, 500, 1000,5000,100_00, 500_00, 1000_00, 1500_00] # mm
	            Ax_2.xticks = (log1p.(Ψticks), string.(cst.Mm_2_kPa .* Ψticks))    
	
	            # ylims!(Ax_2, 0.0, 10.0)
	
	            lines!(Ax_2, Ψ_Log, Ksim, linewidth=Linewidth, color=:red2)
	
	            vlines!(Ax_2, X_Line_θΨ; ymin = 0.0, color=:blue2, linewidth=Linewidth/2.0, linestyle=:dash)
	
	            hlines!(Ax_2, Y_Line_KsMat; xmin=0.0, color=:blue2, linewidth=Linewidth/2.0, linestyle=:dash)
	
	            # scatter!(Ax_2, KΨ_LogΨobs, KΨ_Kobs_1 .* cst.MmS_2_CmH, marker=:hexagon, markersize=25)
	
	            # scatter!(Ax_2, KΨ_LogΨobs, KΨ_Kobs_2 .* cst.MmS_2_CmH, marker=:diamond, markersize=25)
	
	            scatter!(Ax_2, KΨ_Ψobs_Sand, KΨ_Kobs_Sand .* cst.MmS_2_MmH, marker=:diamond, markersize=30, color=:darkgoldenrod3)
	
	            scatter!(Ax_2, KΨ_Ψobs_Silt, KΨ_Kobs_Silt .* cst.MmS_2_MmH  , markersize=30, color=:teal)
	
	            # scatter!(Ax_2, KΨ_Ψobs_Clay, KΨ_Kobs_Clay .* cst.MmS_2_MmH, marker=:diamond, markersize=25)
	
	      # PLOTTING AX_3: DISTRIBUTION
	         Ax_3 = Axis(Fig[1, 1], width=Width, height=Height, xticklabelrotation = π / 4.0, xlabel= L"$R$ [mm]", ylabel="Prob Dens Function", xlabelsize=xlabelSize, ylabelsize=xlabelSize, xgridvisible=xgridvisible, ygridvisible=xgridvisible)
	
	      
	         ylims!(Ax_3, 0.0, 1.2)
	
	         Ax_3.xticks = (log1p.(Ψticks), "10^" .* string.(round.(log10.(cst.Y ./ Ψticks), digits=1)))
	
	         lines!(Ax_3, Ψ_Log, ∂θ∂Ψ_Mat, linewidth=Linewidth, color=:Green)
	         lines!(Ax_3, Ψ_Log, ∂θ∂Ψ_Mac, linewidth=Linewidth, color=:red2)
	
	 Fig
end

# ╔═╡ 57dc88a0-606b-4e0e-8a16-386fd9aedaa0


# ╔═╡ Cell order:
# ╠═777b58c7-180f-41a9-bfa5-edbbc246c939
# ╠═fcfdad85-3e5d-4642-93a7-7843bd2cf50f
# ╠═111ff21d-bc50-4c05-933d-db813261a196
# ╠═119793bf-274a-43b4-9b8c-0aa656adc71d
# ╠═c98dca7b-0612-4c2a-804f-b545acac4b34
# ╠═0e06f203-8f2a-4844-8678-9b2bf9a3f495
# ╠═57dc88a0-606b-4e0e-8a16-386fd9aedaa0
