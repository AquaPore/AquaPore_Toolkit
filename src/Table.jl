# =============================================================
#		MODULE: table
# =============================================================
module table
	import CSV, Tables
	export TABLE_ID

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TABLE_EXTRAPOINTS_K
	# 		Tabular values of the PSD model
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TABLE_ID(NiZ::Int64, path, Path::String)
			println("    ~  $(Path) ~")

			IdSelect = collect(1:1:NiZ)

			Select = fill(1::Int64, NiZ)

			FieldName_String = ["Id", path.option.Select]

			CSV.write(Path, Tables.table([string.(IdSelect[1:NiZ]) Select[1:NiZ]]), writeheader=true, header=FieldName_String, bom=true)
		return nothing
		end  # function:  TABLE_ID


	# =============================================================
	#		MODULE: hydroLab
	# =============================================================
	module hydroLab
		import  ...tool, ...wrc, ...kunsat
		# import DelimitedFiles
		import CSV, Tables
		export θΨK, θΨK_OPTIMISATION

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θΨK(hydro, hydroOther, IdSelect::Vector{Int64}, KₛModel::Vector{Float64}, N_KΨobs::Vector{Int64}, N_θΨobs::Vector{Int64}, NiZ::Int64, Path::String, RockFragment::Vector{Float64})
				println("    ~  $(Path) ~")


				Matrix₁, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(NiZ, hydro)

				Matrix₂, FieldName_String2 = tool.readWrite.STRUCT_2_FIELDNAME(NiZ, hydroOther)

				Header = vcat("Id", FieldName_String, "RockFragment", "KsModel", "N_θΨobs", "N_KΨobs", FieldName_String2 )

				CSV.write(Path, Tables.table([string.(IdSelect) Matrix₁ RockFragment KₛModel N_θΨobs N_KΨobs Matrix₂ ]), writeheader=true, header=Header, bom=true)
			return nothing
			end  # function:  θΨK
		#---------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK_OPTIMISATION
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θΨK_OPTIMISATION(hydro, hydroOther, NiZ, Path)

				Matrix₁, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(NiZ, hydro)

				Matrix₂, FieldName_String2 = tool.readWrite.STRUCT_2_FIELDNAME(NiZ, hydroOther)

				Header = vcat( FieldName_String, FieldName_String2 )

				CSV.write(Path, Tables.table([ Matrix₁ Matrix₂ ]), writeheader=true, header=Header, bom=true)
			return nothing
			end  # function:  θΨK
		#---------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : TABLE_θΨK
		# 		Tabular values of the PSD model
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function TABLE_EXTRAPOINTS_θΨ(optionₘ, hydro, IdSelect, NiZ::Int64, Path, Ψ_Table; Orientation="Horizontal")
				println("    ~  $(Path) ~")

				N_Ψ = Int64(length(Ψ_Table))

				if Orientation == "Horizontal" # <>=<>=<>=<>=<>=<>
					# Writting the Header
						FieldName_String = fill(""::String, N_Ψ)
						for i =1:N_Ψ
							FieldName_String[i] = string(Int64(Ψ_Table[i]) ) * "mm"
						end
						pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
							
					# Computing θ at required θ
						θ₂ = fill(0.0::Float64, (NiZ, N_Ψ))

						for iZ=1:NiZ
							for iΨ =1:N_Ψ
								Ψ₂ = Ψ_Table[iΨ]
								θ₂[iZ, iΨ] = wrc.Ψ_2_θ(optionₘ, Ψ₂, iZ, hydro)
							end # iΨ
						end # iZ

						CSV.write(Path, Tables.table([Int64.(IdSelect) θ₂]), writeheader=true, header=FieldName_String, bom=true)

				elseif Orientation == "Vertical" # <>=<>=<>=<>=<>=<>
					FieldName_String = ["Id","Ψ[mm]","Theta[0-1]"]
					N = N_Ψ * NiZ
					Id₂ = fill(0::Int64, N)
					Ψ₂  = fill(0.0::Float64, N)
					θ₂  = fill(0.0::Float64, N)
					iCount = 1

					for iZ=1:NiZ
						for iΨ =1:N_Ψ
							Id₂[iCount] = IdSelect[iZ]
							Ψ₂[iCount] = Ψ_Table[iΨ]
							θ₂[iCount] = wrc.Ψ_2_θ(optionₘ,Ψ₂[iCount], iZ, hydro)
	
							iCount+=1
						end # iΨ
					end # iZ

					CSV.write(Path, Tables.table([string.(Id₂[1:N]) Ψ₂[1:N] θ₂[1:N]]), writeheader=true, header=FieldName_String, bom=true)
				else
					error("SoilWaterToolBox Error in TABLE_EXTRAPOINTS_θΨ $Orientation not understood")
				end
		return nothing
		end  # function:  θΨ

		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : TABLE_EXTRAPOINTS_K
		# 		Tabular values of the PSD model
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function TABLE_EXTRAPOINTS_Kθ(optionₘ, hydroParam, IdSelect, K_Table, NiZ::Int64, Path::String)
				println("    ~  $(Path) ~")

				N_K = Int64(length(K_Table))

			# Writting the Header
				FieldName_String =["Id", "H[mm]" ,"Kunsat[mm_s]"]
							
			# Computing K at required Ψ
				N = N_K * NiZ
				Id₂     = fill(0::Int64, N)
				Ψ₂      = fill(0.0::Float64, N)
				Kunsat₂ = fill(0.0::Float64, N)
				iCount  = 1
				hydroParam₂ = deepcopy(hydroParam)

				for iZ=1:NiZ
					for iK =1:N_K
						Id₂[iCount] = IdSelect[iZ]
						Ψ₂[iCount] = K_Table[iK]
						Kunsat₂[iCount] = kunsat.KUNSAT_θΨSe(optionₘ, Ψ₂[iCount], iZ, hydroParam₂)

						iCount += 1
					end # iΨ
				end # iZ
				CSV.write(Path, Tables.table([string.(Id₂[1:N]) Ψ₂[1:N] Kunsat₂[1:N]]), writeheader=true, header=FieldName_String, bom=true)
		return nothing
		end  # function:  θΨ
			
	end  # module hydro
	# ............................................................


	# =============================================================
	#		module: ksmodel
	# =============================================================
	module ksmodel
		import ...tool, ...cst
		# import DelimitedFiles
		import CSV, Tables
		export KSMODEL_τ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KSMODEL
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function KSMODEL(hydro, IdSelect, KₛModel, Path)
				println("    ~  $Path ~")

				# Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(2,  ksmodelτ)
				Header = ["IdSelect", "θs", "θsMacMat","θr","Ψm[mm]","LnΨm[mm]","σ","σMac","ΨmMac","ΔθsMacMat-θr","Δθs-ΔθsMacMat","Ks[mm h⁻¹]", "KsModel[mm h⁻¹]","LnKs[mm h⁻¹]","LnKsModel[mm h⁻¹]","ΔKs-KsModel[mm h⁻¹]","ΔLnKs-LnKsModel[mm h⁻¹]"]

				KsObs = hydro.Ks .* cst.MmS_2_MmH
				kₛ_Model₂ = KₛModel .* cst.MmS_2_MmH

				LnKsObs = log1p.(KsObs)
				Lnkₛ_Model = log1p.(kₛ_Model₂)

				ΔKsKsModel = KsObs .- kₛ_Model₂
				ΔLnKsKsModel = LnKsObs .- Lnkₛ_Model

				CSV.write(Path, Tables.table([string.(Int64.(IdSelect)) hydro.θs hydro.θsMacMat hydro.θr hydro.Ψm log.(hydro.Ψm) hydro.σ hydro.σMac hydro.ΨmMac hydro.θsMacMat.-hydro.θr hydro.θs.-hydro.θsMacMat KsObs kₛ_Model₂ LnKsObs Lnkₛ_Model ΔKsKsModel ΔLnKsKsModel ]), writeheader=true, header=Header, bom=true)	
			return nothing
			end  # function: KSMODEL
		# ------------------------------------------------------------------


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : KSMODEL_τ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function KSMODEL_τ(ksmodelτ, Ngroup, Path)
				println("    ~  $Path ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(Ngroup, ksmodelτ)
				
				X, Y = size(Matrix)
				CSV.write(Path, Tables.table(round.(Matrix, digits=5)), writeheader=true, header=FieldName_String, bom=true)	
			return nothing
			end  # function: KSMODEL_τ
			# ------------------------------------------------------------------
		
	end  # module: ksmodel
	# ............................................................


	# =============================================================
	#		MODULE: psd
	# =============================================================
	module psd
		import...tool, ...wrc, ...cst
		import DelimitedFiles
		import CSV, Tables
		export PSD, θΨK_PSD

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PSD
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PSD(IdSelect, NiZ::Int64, paramPsd, Path)
				println("    ~  $Path ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(NiZ,  paramPsd)
				
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

				# open(Path, "w") do io
				# 	DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
				# 	DelimitedFiles.writedlm(io, [string.(Int64.(IdSelect)) round.(Matrix,digits=5)], ",")
				# end
				CSV.write(Path, Tables.table([string.(Int64.(IdSelect)) round.(Matrix,digits=5)]), writeheader=true, header=FieldName_String, bom=true)
			return nothing
			end


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK_PSD
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θΨK_PSD(hydroPsd, IdSelect, KunsatModel_Psd, NiZ::Int64, Path)
				println("    ~  $Path ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(NiZ, hydroPsd)

				Matrix = hcat(Matrix, KunsatModel_Psd)

				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
				push!(FieldName_String, "Table_KΨ")

				Matrix =  round.(Matrix, digits=5)
				# open(Path, "w") do io
				# 	DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
				# 	DelimitedFiles.writedlm(io, [string.(IdSelect) Matrix], ",")
				# end
				CSV.write(Path, Tables.table([string.(IdSelect) Matrix]), writeheader=true, header=FieldName_String, bom=true)
			return nothing
			end  # function:  θΨK_PSD



		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PSD_θΨ_θ
		# 		Tabular values of the PSD model
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PSD_θΨ_θ(IdSelect, hydroPsd, NiZ::Int64, option, param, Path::String)
				println("    ~  $Path ~")

				N_Ψ = Int64(length(param.psd.Ψ_Table))

				# Writting the Header
					FieldName_String = fill(""::String, (N_Ψ))

					for i =1:N_Ψ
						FieldName_String[i] = string(param.psd.Ψ_Table[i] * cst.Mm_2_Cm) * "cm"
					end
					pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
				
				# Computing θ at required θ
					θ = fill(0.0::Float64, (NiZ, N_Ψ))
					for iZ=1:NiZ
						for iΨ =1:N_Ψ
							Ψ = param.psd.Ψ_Table[iΨ]
							θ[iZ, iΨ] = wrc.Ψ_2_θ(option.psd, Ψ, iZ, hydroPsd)
						end # iΨ
					end # iZ

				# Concatenating the 2 matrices
					θ = hcat(IdSelect, θ)
					θ = round.(θ, digits=5)

				# Writting the table
					open(Path, "w") do io
						DelimitedFiles.writedlm(io, [FieldName_String] , ",",) # Header
						for i = 1:length(IdSelect)
							DelimitedFiles.writedlm(io, [θ[i, 1:N_Ψ+1]], ",")
						end # i
					end # Path
			return nothing
			end  # function:  θΨK_PSD

	end  # module psd
	# ............................................................
	
	

	# =============================================================
	#		MODULE: infilt
	# =============================================================
	module infilt
		import ...tool
		import CSV, Tables
		export HYDRO_INFILT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYDRO_INFILT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYDRO_INFILT(hydroInfilt, IdSelect, NiZ::Int64, Path, Soilname)
				println("    ~  $Path ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(NiZ, hydroInfilt)

				# Matrix = hcat(Matrix, KunsatModel_Infilt)\

				pushfirst!(FieldName_String, "Soilname") # Write the "Soilname" at the very begenning
				# push!(FieldName_String, string("Kunsat_θΨ"))
				
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
				# push!(FieldName_String, string("Kunsat_θΨ"))

				Matrix =  round.(Matrix, digits=10)
				CSV.write(Path, Tables.table( [string.(IdSelect) Soilname Matrix]), writeheader=true, header=FieldName_String, bom=true)
			return nothing
			end  # function: HYDRO_INFILT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : infilt
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function INFILT(IdSelect, infiltOutput, NiZ::Int64, Path, Soilname)
				println("    ~  $Path ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(NiZ::Int64, infiltOutput)

				pushfirst!(FieldName_String, "Soilname") # Write the "Soilname" at the very begenning				
				
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

				Matrix =  round.(Matrix, digits=5)
				CSV.write(Path, Tables.table( [string.(IdSelect) Soilname Matrix]), writeheader=true, header=FieldName_String, bom=true)
			return nothing
			end  # function: HYDRO_INFILT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : infilt
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function INFILT_SEINI(IdSelect, Infilt_SeIni, NiZ::Int64, Path, Soilname, Time_2_Infilt)
				println("    ~  $Path ~")

				FieldName_String = "SeIni=" .* string.(Infilt_SeIni)

				pushfirst!(FieldName_String, "Soilname") # Write the "Soilname" at the very begenning					

				pushfirst!(FieldName_String, string("Id/Seconds")) # Write the "Id" at the very begenning

				Matrix =  round.(Time_2_Infilt, digits=1)
				CSV.write(Path, Tables.table( [string.(IdSelect) Soilname Matrix]), writeheader=true, header=FieldName_String, bom=true)
			return nothing
			end  # function: HYDRO_INFILT
		
	end  # module: infilt

	# .........................................................................................


	# =============================================================
	#		module: other
	# =============================================================
		module convert
			import CSV, Tables
			export CONVERT_θΨ_2D_2_1D
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : CONVERT_θΨ_2D_2_1D
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function CONVERT_θΨ_2D_2_1D(IdSelect, NiZ, N_θΨobs, Path, θ_θΨobs, Ψ_θΨobs)
					println("			~  $(Path) ~")

					Header = ["Id","H[mm]","Theta[0-1]"]

					Ψ_1D=[]; θ_1D=[]; Id_Repeat=[]

					for iZ = 1:NiZ
						for iΨ = 1:N_θΨobs[iZ]
							append!(Id_Repeat, IdSelect[iZ])
							append!(Ψ_1D, Ψ_θΨobs[iZ, iΨ])
							append!(θ_1D, θ_θΨobs[iZ, iΨ])
						end # for: iΨ = 1:Ψ_θΨobs
					end # for: iZ = 1:NiZ

				# open(Path, "w") do io
				# 	DelimitedFiles.writedlm(io,[Header] , ",",) # Header
				# 	DelimitedFiles.writedlm(io, [Id_Repeat Ψ_1D θ_1D] , ",")
				# end # open
				CSV.write(Path, Tables.table([Id_Repeat Ψ_1D θ_1D]), writeheader=true, header=Header, bom=true)
			return nothing
			end  # function: CONVERT_θΨ_2D_2_1D
			# ------------------------------------------------------------------


			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : CONVERT_KΨ_2D_2_1D
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function CONVERT_KΨ_2D_2_1D(IdSelect, NiZ, N_KΨobs, Path, K_KΨobs, Ψ_KΨobs)
					println("			~  $(Path) ~")

					Header = ["Id","H[mm]","Kunsat[mm_s]"]

					Ψ_1D=[]; θ_1D=[]; Id_Repeat=[]

					for iZ = 1:NiZ
						for iΨ = 1:N_KΨobs[iZ]
							append!(Id_Repeat, IdSelect[iZ])
							append!(Ψ_1D, Ψ_KΨobs[iZ, iΨ])
							append!(θ_1D, K_KΨobs[iZ, iΨ])
						end # for: iΨ = 1:Ψ_θΨobs
					end # for: iZ = 1:NiZ

				# open(Path, "w") do io
				# 	DelimitedFiles.writedlm(io,[Header] , ",",) # Header
				# 	DelimitedFiles.writedlm(io, [Id_Repeat Ψ_1D θ_1D] , ",")
				# end # open
				CSV.write(Path, Tables.table( [Id_Repeat Ψ_1D θ_1D]), writeheader=true, header=Header, bom=true)
			return nothing
			end  # function: CONVERT_KΨ_2D_2_1D
		# ------------------------------------------------------------------
				
		end  # module: convert
		# ............................................................

end  # module table
# ............................................................