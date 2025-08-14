library(profvis)

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}


#full testing

a <- Single_Plot(Data = Intelligence_Sum_Stats, Verbose = TRUE)
b <- Regional_Plot(Data = Intelligence_Sum_Stats, Auto_LD = F, Region_Window = 1e5 , Chromosomes = c(1:5), Verbose = T)
c <- Miami_Plot(Top_Data = Intelligence_Sum_Stats, Bottom_Data = Intelligence_Sum_Stats)


Regional_Plots <- Regional_Plot(Data = Intelligence_Sum_Stats, Chromosomes = 1, Verbose = T)

 b <- Regional_Plot(Data = Intelligence_Sum_Stats,  Chromosomes = 1, Auto_LD = FALSE, Region_Window = 1e5, Verbose = T )


profvis({
 a <-  Single_Plot(Data = Intelligence_Sum_Stats)
  #  "C:/Users/callumon/Miami_Package_R/MiamiR/Intelligence_Sum_Stats_Mini.txt"
})



a <- run_with_counter(
  Regional_Plot,
  args = list(
    Data = Intelligence_Sum_Stats,
    .dots = list()  # <- ensures .dots is defined
  ),
  default_val = NULL
)



a <- run_with_counter(Regional_Plot, args = list(Data = Intelligence_Sum_Stats))


a <- run_with_counter(Miami_Plot, args = list(Top_Data = Intelligence_Sum_Stats,
                                              Bottom_Data = Intelligence_Sum_Stats))



profvis({
  Single_Plot(Data = Intelligence_Sum_Stats)
  #  "C:/Users/callumon/Miami_Package_R/MiamiR/Intelligence_Sum_Stats_Mini.txt"
})


profvis({
  Miami_Plot(Intelligence_Sum_Stats)
  #  "C:/Users/callumon/Miami_Package_R/MiamiR/Intelligence_Sum_Stats_Mini.txt"
})


a <- Single_Plot(Data = Intelligence_Sum_Stats)











