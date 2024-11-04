
module evapoIran
using CSV, Tables, DataFrames
using CairoMakie, ColorSchemes
using BlackBoxOptim
include(raw"D:\MAIN\MODELS\AquaPore_Toolkit\src\Stats.jl")

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : EVAPOTRANSPIRATION_IRAN
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function EVAPOTRANSPIRATION_IRAN(;Model, MaxFuncEvals)

		println("==============================================")
		println("=== Model=$Model ====")

		# Parameters
         P₁_Max = 500.0; P₁_Min = 10.0
         P₂_Max = 10.0;  P₂_Min = 1.0
         P₃_Max = 500.0;    P₃_Min = 0.1
         P₄_Max = 500.0;    P₄_Min = 0.1
         P₅_Max = 1.;  P₅_Min = 0.0

			NparamOpt=1
			if Model == 1
				NparamOpt = 1

			elseif Model == 2
				NparamOpt = 2

			elseif Model == 3
				NparamOpt = 3

			elseif Model == 4
				NparamOpt = 3

			elseif Model == 5
				NparamOpt = 4
			end

			X= zeros(Float64,5)
			X[1]=P₁_Min
			X[2]=P₂_Min
			X[3]=P₃_Min
			X[4]=P₄_Min
			X[5]=P₅_Min

		# Reading data from csv
			Path_Input = raw"D:\MAIN\MODELS\AquaPore_Toolkit\data\INPUT\Data_Evapotranspiration\EvapotranspitationIran.csv"

         Data      = CSV.read(Path_Input, DataFrame;  missingstring=["NA", "NAN", ""], ignoreemptyrows=true)

         Se        = Vector(Data[!,"Se"])
         Clay      = Vector(Data[!,"Clay"])
         Sand      = Vector(Data[!,"Sand"])
         Om        = Vector(Data[!,"OM"])
         Evapo_Obs = Vector(Data[!,"EvapoObs"])
		
		# Initialising
			N = length(Se)
			Evapo_Sim = zeros(Float64, N)

			println("Optimising NparamOpt=$NparamOpt")

			SearchRange = SEARCHRANGE(NparamOpt, P₁_Max, P₁_Min, P₂_Max, P₂_Min, P₃_Max, P₃_Min, P₄_Max, P₄_Min, P₅_Max, P₅_Min)
			println("\n === Feasible range ===")
			println(SearchRange)
			println("\n")

		# Optimisation
			Optimization = BlackBoxOptim.bboptimize(X ->  evapoIran.ObjectiveFunction(Clay, Evapo_Obs, Evapo_Sim, Model, N, NparamOpt, Om, Sand, Se, X); SearchRange=SearchRange, NumDimensions=NparamOpt, TraceMode=:silent)

			# , MaxFuncEvals=MaxFuncEvals

         X            = BlackBoxOptim.best_candidate(Optimization)
         Of           = evapoIran.ObjectiveFunction(Clay, Evapo_Obs, Evapo_Sim, Model, N, NparamOpt, Om, Sand, Se, X)
         Evapo_Sim = evapoIran.EVAPOTRANSPIRATION_MODEL(Clay, Evapo_Sim, Model, N, Om, Sand, Se, X)

		# Output
			println("Model=$Model, NumDimensions=$NparamOpt, Of=$Of")
			println("X parameters =")
			@show X
		
			Rmse = stats.RMSE(Evapo_Obs[1:N], Evapo_Sim[1:N])
			Nse = stats.NSE(Evapo_Obs[1:N], Evapo_Sim[1:N])

			println("\n RMSE = $Rmse")
			println("NSE = $Nse")

		# Plotting
			evapoIran.PLOT_EVAPO(Evapo_Obs[1:N], Evapo_Sim[1:N], Model)

	return nothing
	end  # function: EVAPOTRANSPI	RATION_IRAN
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : EVAPOTRANSPIRATION_MODEL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function EVAPOTRANSPIRATION_MODEL(Clay, Evapo_Sim, Model, N, Om, Sand, Se, X)

			@inbounds for i =1:N
				if Model == 1
					Evapo_Sim[i] = X[1] * sin(Se[i] * π * 0.5)

				elseif Model == 2
					Evapo_Sim[i] = X[1] * sin(Se[i] * π * 0.5) ^ X[2]
			
				elseif Model == 3
					Evapo_Sim[i] = X[1] * sin(Se[i] * π * 0.5) ^ max(X[2] - X[3] * Clay[i], 0.5)

				elseif Model == 4
					Evapo_Sim[i] = X[1] * sin(Se[i] * π * 0.5) ^ (X[2] + X[3] * Sand[i])

				elseif Model == 5
					Evapo_Sim[i] = X[1] * sin(Se[i] * π * 0.5) ^ max(X[2] - X[3] * Clay[i] + X[4] * Sand[i], 0.5)
				end
			end		
		return Evapo_Sim
		end  # function: EVAPOTRANSPIRATION_MODEL
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ObjectiveFunction
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ObjectiveFunction(Clay, Evapo_Obs, Evapo_Sim, Model, N, NparamOpt, Om, Sand, Se, X)

			Evapo_Sim = evapoIran.EVAPOTRANSPIRATION_MODEL(Clay, Evapo_Sim, Model, N, Om, Sand, Se, X)

		return 1.0 - abs(stats.NSE_WILMOT(Evapo_Obs, Evapo_Sim))
		end  # function: ObjectiveFunction
	# ------------------------------------------------------------------


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SEARCHRANGE
	#		Required by BlackBoxOptim
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SEARCHRANGE(NparamOpt, P₁_Max, P₁_Min, P₂_Max, P₂_Min, P₃_Max, P₃_Min, P₄_Max, P₄_Min, P₅_Max, P₅_Min)
			ParamOpt_Min₂ = [P₁_Min, P₂_Min, P₃_Min, P₄_Min, P₅_Min]
			ParamOpt_Max₂ = [P₁_Max, P₂_Max, P₃_Max, P₄_Max, P₅_Max]
		return SearchRange = (collect(zip(Float64.(ParamOpt_Min₂[1:NparamOpt]), Float64.(ParamOpt_Max₂[1:NparamOpt]))))
		end  # function: SEARCHRANGE
	#..................................................................
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT_EVAPO
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Plotting parameters
			ColourOption_No    = 1
			Linewidth          = 2
			height             = 200
			labelsize          = 15
			Markersize         = 8
			textcolor          = :blue
			textsize           = 20
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
			xticklabelrotation = 0
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

		function PLOT_EVAPO(Evapo_Obs, Evapo_Sim, Model)

			CairoMakie.activate!(type="svg", pt_per_unit=1)
			Fig =  Figure(figure_padding = 10; fonts = ( ; regular="CMU Serif")) 

			Axis_ObsSim = Axis(Fig[1,1], xlabel= L"$Evapot Obs$ [mm]", ylabel=L"$Evapot Sim$  [mm]", title= "Model" * string(Model), titlecolor=titlecolor, xticklabelrotation=xticklabelrotation, ylabelsize=ylabelsize, xlabelsize=xlabelSize, xticksize=xticksize, yticksize=yticksize, width=width, height=height,   titlesize=titlesize,  xgridvisible=xgridvisible, ygridvisible=ygridvisible, xminorticksvisible=xminorticksvisible, yminorticksvisible=yminorticksvisible, xtickwidth=xtickwidt, ytickwidth=ytickwidt, xtickalign=xtickalign, ytickalign=ytickalign, xticksmirrored=xticksmirrored, yticksmirrored=yticksmirrored,  xtrimspine=xtrimspine,  ytrimspine=ytrimspine, xgridstyle=xgridstyle, ygridstyle=ygridstyle, yminorticks=IntervalsBetween(5), xlabelpadding=xlabelpadding, ylabelpadding=ylabelpadding, xminortickalign=xminortickalign, yminortickalign=yminortickalign)

			scatter!(Axis_ObsSim, Evapo_Obs, Evapo_Sim, markersize=Markersize, marker = '●', strokewidth=1, strokecolor=:deepskyblue3, label="")

				Evapo_Max = max(maximum(Evapo_Obs), maximum(Evapo_Sim))
				Line = range(0.0, stop=Evapo_Max, length=10) 
				lines!(Axis_ObsSim, Line, Line, color=:grey, linestyle=:dash, linewidth=2)

			# General
				resize_to_layout!(Fig)
				trim!(Fig.layout)
				colgap!(Fig.layout, 10)
				rowgap!(Fig.layout, 10)

				Path =raw"D:\MAIN\MODELS\AquaPore_Toolkit\data\OUTPUT\Evapotranspiration\Evapotranspiration_" * string(Model) * ".svg"
				save(Path, Fig)
				display(Fig)	
				
			return nothing
			end  # function: PLOT_EVAPO
		# ------------------------------------------------------------------
end

evapoIran.EVAPOTRANSPIRATION_IRAN(Model=1, MaxFuncEvals=10000)
evapoIran.EVAPOTRANSPIRATION_IRAN(Model=2, MaxFuncEvals=10000)
evapoIran.EVAPOTRANSPIRATION_IRAN(Model=3, MaxFuncEvals=10000)
evapoIran.EVAPOTRANSPIRATION_IRAN(Model=4, MaxFuncEvals=10000)
evapoIran.EVAPOTRANSPIRATION_IRAN(Model=5, MaxFuncEvals=10000)