
profvis::profvis({

  Manhattan_Plot <- Single_Plot(Data = Intelligence_Sum_Stats, Verbose = T)

})

profvis::profvis({

  Regional_Plot <- Regional_Plot(Data = Household_Income_Sum_Stats, Chromosomes = 1, Verbose = T)

})

profvis::profvis({

  Miami_Plot <- Miami_Plot(Top_Data = Household_Income_Sum_Stats, Bottom_Data = Intelligence_Sum_Stats, Verbose = T)

})


