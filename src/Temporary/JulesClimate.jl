# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : JULES_CLIMATE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import CSV, Tables

function JULES_CLIMATE()

   Path₀ = "D:\\DATAraw\\JULES_DATA\\NO_STONES"

   cd(Path₀)

   ClimateInventory = readdir()

   println(ClimateInventory)
   i=0

   ClimatePath_Output = joinpath(raw"D:\DATAraw\JULES_DATA\DATE\DATE_SUMMARY.csv" )

   try
      rm(ClimatePath_Output)
   catch
      nothing
   end

   for  iClimateInventory in ClimateInventory
      i+=1
      ClimatePath_Input = joinpath(Path₀, iClimateInventory, iClimateInventory * "_Dates.csv" ) 

      ClimatePath_Output = joinpath(raw"D:\DATAraw\JULES_DATA\DATE", iClimateInventory * "_Date.csv" )


      # cp(ClimatePath_Input, ClimatePath_Output)

      Data = CSV.File(ClimatePath_Input, header=true)

      Year_Obs_Start  = convert(Vector{Float64}, Tables.getcolumn(Data, :Year_Obs_Start))
      Month_Obs_Start = convert(Vector{Float64}, Tables.getcolumn(Data, :Month_Obs_Start))
      Day_Obs_Start   = convert(Vector{Float64}, Tables.getcolumn(Data, :Day_Obs_Start))

      Year_Obs_End    = convert(Vector{Float64}, Tables.getcolumn(Data, :Year_Obs_End))
      Month_Obs_End   = convert(Vector{Float64}, Tables.getcolumn(Data, :Month_Obs_End))
      Day_Obs_End     = convert(Vector{Float64}, Tables.getcolumn(Data, :Day_Obs_End))

      Year_Sim_Start  = convert(Vector{Float64}, Tables.getcolumn(Data, :Year_Sim_Start))
      Month_Sim_Start = convert(Vector{Float64}, Tables.getcolumn(Data, :Month_Sim_Start))
      Day_Sim_Start   = convert(Vector{Float64}, Tables.getcolumn(Data, :Day_Sim_Start))

      Year_Sim_Start  = convert(Vector{Float64}, Tables.getcolumn(Data, :Year_Sim_Start))
      Month_Sim_Start = convert(Vector{Float64}, Tables.getcolumn(Data, :Month_Sim_Start))
      Day_Sim_Start   = convert(Vector{Float64}, Tables.getcolumn(Data, :Day_Sim_Start))

      Year_Sim_End    = convert(Vector{Float64}, Tables.getcolumn(Data, :Year_Sim_End))
      Month_Sim_End   = convert(Vector{Float64}, Tables.getcolumn(Data, :Month_Sim_End))
      Day_Sim_End     = convert(Vector{Float64}, Tables.getcolumn(Data, :Day_Sim_End))

      Year_Sim_Start  = Year_Obs_Start .+ 1
      Month_Sim_Start = Month_Obs_Start
      Day_Sim_Start   = Day_Obs_End

      Hour_Start=[9]
      Hour_Start_Sim=[9]
      Hour_End =[9]

      NEW_FILE(iClimateInventory, Year_Obs_Start, Year_Sim_Start, Year_Obs_End ,Month_Obs_Start, Month_Sim_Start ,Month_Obs_End, Day_Obs_Start ,Day_Sim_Start ,Day_Obs_End ,Hour_Start,Hour_Start_Sim ,Hour_End, i)
   end
   
return nothing
end  # function: JULES_CLIMATE
# ------------------------------------------------------------------


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : NEW_FILE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function NEW_FILE(iClimateInventory, Year_Start, Year_Start_Sim, Year_End ,Month_Start, Month_Start_Sim ,Month_End, Day_Start ,Day_Start_Sim ,Day_End ,Hour_Start, Hour_Start_Sim ,Hour_End, i)

   ClimatePath_Output = joinpath(raw"D:\DATAraw\JULES_DATA\DATE\DATE_SUMMARY.csv" )
      
   Header = [ "Station", "Year_Start" ,"Year_Start_Sim", "Year_End" , "Month_Start" ,"Month_Start_Sim" ,"Month_End" ,"Day_Start", "Day_Start_Sim" ,"Day_End" ,"Hour_Start" ,"Hour_Start_Sim", "Hour_End"]

   CSV.write(ClimatePath_Output, Tables.table([iClimateInventory Year_Start Year_Start_Sim Year_End Month_Start Month_Start_Sim Month_End Day_Start Day_Start_Sim Day_End Hour_Start Hour_Start_Sim Hour_End]), bom=true, append=true,   writeheader = true,header=Header, delim = ",")

      
   return nothing
end  # function: NEW_FILE
# ------------------------------------------------------------------

JULES_CLIMATE()


