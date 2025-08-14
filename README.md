---
title: "README"
output: html_document
---

```{r, include = FALSE}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

```

Welcome to the MiamiR software package overview/tutorial - key for all your GWAS and other genetic data visualisation needs!

Firstly, the recommended installation for this package is as follows: remotes::install_github("CRON-Project/MiamiR") - Alternative method: BiocManager::install("CRON-Project/MiamiR")

If you have any issues please contact me at: callum.oneill@rocketmail.com with subject heading indicating issue and please raise issues through GitHub!

Let's load in the MiamiR package (lib.loc may need to be specified), which is complete with extensive documentation and functionality. A web app version is also available at: TBC

```{r setup, message=TRUE}

library(MiamiR)

```

Here's how to create a single Manhattan plot of the included 'Intelligence_Sum_Stats' summary statistics data set using the Single_Plot() function from MiamiR.

This function will return a plot object which we'll manually save to a .jpg file in the vignettes/Plots directory of the package for future use. Initially, all automated, default settings will be used.

```{r, results="hide", warning = FALSE, message = TRUE}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

Manhattan_Plot <- Single_Plot(Data = Intelligence_Sum_Stats)

ggplot2::ggsave("Intelligence_Plot.jpg", plot = Manhattan_Plot, width = 30,
                 height = 15, units = "in", dpi = 100)

```
Let's inspect our saved plot.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=300, out.width="100%"}

knitr::include_graphics("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots/Intelligence_Plot.jpg")

```

The Single_Plot() function from MiamiR also allows for a range of customisation. Here is a limited example of such functionality.

Note: running ?Single_Plot() or ? for any function in MiamiR will provide documentation on all functionality available and information on any in-house data sets used.

Also real-time progress bars and outputs are available for all functions. Please specify Verbose = TRUE, for extended versions of real-time checkpoints, included mainly for debugging and reporting purposes.

```{r, results="hide", warning = FALSE, message = TRUE}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

Manhattan_Plot <- Single_Plot(Data = Intelligence_Sum_Stats, X_Axis_Title = "Region", Title = "Custom Plot", Label_Height = 15, Anchor_Label = "centre", 
                              Label_Angle = 90, Chromosome_Colours = c("red", "orange"), Colour_Of_Index = "blue", Diamond_Index = T, Diamond_Index_Size = 10)

ggplot2::ggsave("Intelligence_Plot_Custom.jpg", plot = Manhattan_Plot, width = 30,
                 height = 15, units = "in", dpi = 100)

```
Let's inspect our saved plot.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=300, out.width="100%"}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

knitr::include_graphics("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots/Intelligence_Plot_Custom.jpg")

```

Files can also be fed in for all MiamiR functions using variables containing a file path or a raw file path, in addition to multiple files or environment objects, whilst allowing for identical functionality.

Here, raw file paths for the above summary statistics are used alongside another included data set 'Household_Income_Summary_Statistics', with extra customisation demonstrated as well.

Notice how automated, unspecified defaults also subtly adjust for different inputs.

```{r, results="hide", warning = FALSE, message = TRUE}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

Sum_Stats <- c("C:/Users/callumon/Miami_Package_R/MiamiR/inst/extdata/Example_Data_Raw/Intelligence_Sum_Stats.txt",
               "C:/Users/callumon/Miami_Package_R/MiamiR/inst/extdata/Example_Data_Raw/Household_Income_Sum_Stats.txt")

Manhattan_Plot <- Single_Plot(Data = Sum_Stats, Label_Angle = 60, Label_Colour = "grey", Sig_Line_Type = "solid", Condense_Scale = F )

#Save names as looping to display later

output_files <- c()
for (plot_name in names(Manhattan_Plot)) {
  file_path <- paste0(plot_name, ".jpg")
  ggplot2::ggsave(
    filename = file_path,
    plot = Manhattan_Plot[[plot_name]],
    width = 30, height = 15, units = "in", dpi = 100
  )
  output_files <- c(output_files, file_path)
}

```
Let's inspect our saved plots.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")
knitr::include_graphics(output_files)

```

There is also a Regional_Plot() function built around the above functionality, which has unique features in addition to seamlessly inheriting the Single_Plot() arguments to produce annotated regional plots. 

Here is how to run the function, which automatically ascertains regions of interest, with this example focusing on key regions in chromosomes 10 and 22 of the intelligence GWAS summary statistics.

We will also save all generated objects as correctly sized images with guide measures also returned by the function and supplied as part of the plot objects. 

Although defaults are suitable it is very important to specify the correct genomic build (in this case HG19) for accurate recombination and gene annotation tracks/overlays, in addition to queried LD, where needed.

```{r,  results="hide",  warning = FALSE, message = TRUE}

  setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

  Regional_Plots <- Regional_Plot(Data = Intelligence_Sum_Stats,  Genome_Build = "grch37", Chromosomes = c(10,22))
  
  for (plot_name in names(Regional_Plots)) {
  plot <- Regional_Plots[[plot_name]]
  file_safe <- gsub("[^A-Za-z0-9]", "_", plot_name)
  filename <- paste0(file_safe, ".jpg")
  height_to_use <- attr(plot, "dynamic_height")
  if (is.null(height_to_use)) height_to_use <- 30

  setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")
  
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = 30,
    height = height_to_use,  
    units = "in",
    dpi = 100,
    limitsize = FALSE
  )
  
}

```
Let's inspect our regional plots for the aforementioned chromosomes.
 
```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

  setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")
  
  img_dir <- "C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots"
  
  # List files that end with 'Regional_Plot.jpg'
  img_paths <- list.files(
    path = img_dir,
    pattern = "Regional_Plot\\.jpg$",  # exact match
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Sort by modification time (oldest first)
  img_paths <- img_paths[order(file.info(img_paths)$mtime)]
  
  # Display all matching plots
  knitr::include_graphics(img_paths)

```

One key customisation of the Regional_Plot() function is the ability to query for LD values around the lead SNP (with settings available for build, population etc.) and colour points in the top panel accordingly. 

For simplicity and speed let's focus on the single chromosome 22 peak. Please bear in mind that querying an API is slow and may not always contain data pertaining to your particular set of SNPs. 

It is often better to produce your own LD matrix using tools such as PLINK1.9/2, the results of which can also be fed in to this package (FUTURE/TBC).

```{r, results="hide", warning = FALSE, message = TRUE}

  setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

  Regional_Plots <- Regional_Plot(Data = Sum_Stats[1], Auto_LD = TRUE, Chromosomes = c(22), Genome_Build = "grch37", 
                                  Gene_Biotype = c("coding", "snorna"),  Gene_Biotype_Colour = c("red", "green" ),
                                  Sense_Arrow_Colour = "darkorange", Gene_Structure_Size_Exons = 20, Gene_Legend_Location = "Top Left",
                                  Gene_Legend_Title = "Categories", Y_Axis_Title = "LOG10Pvalue", Diamond_Index_Size = 20)

  setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

  for (plot_name in names(Regional_Plots)) {
  plot <- Regional_Plots[[plot_name]]
  file_safe <- gsub("[^A-Za-z0-9]", "_", plot_name)
  filename <- paste0(file_safe, "_LD", ".jpg")
  height_to_use <- attr(plot, "dynamic_height")
  if (is.null(height_to_use)) height_to_use <- 30

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = 30,
    height = height_to_use,   #nice divide by exon size to get ratio
    units = "in",
    dpi = 100,
    limitsize = FALSE
  )
  
}

```
Let's have a look at our heavily customised and LD annotated regional plot.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

  setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")
  
  img_dir <- "C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots"
  
  # List files that end with 'Regional_Plot.jpg'
  img_paths <- list.files(
    path = img_dir,
    pattern = "Regional_Plot_LD\\.jpg$",  # exact match
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Sort by modification time (oldest first)
  img_paths <- img_paths[order(file.info(img_paths)$mtime)]
  
  # Display all matching plots
  knitr::include_graphics(img_paths)
  
```

Also, let's now pass both sets of summary statistics we used in the earlier example to show how regional plots can be constructed on distinct data sets simultaneously.

Both Intelligence_Sum_Stats and Household_Income_Sum_Stats seem to have a single significant association on chromosome 9 - let's see if they are in similar genomic locations...

For single vector arguments like Point_Colour, a list can be allocated in the order of customisation of data required. 

Point_Colour = c("blue", "red") will allocate blue to the points in Intelligence_Sum_Stats and red to the points in Household_Income_Sum_Stats. Any singularly specified arguments will apply globally.
 
Arguments which require a list structure, like Gene_Biotype can have two separate listed choices passed by allotting e.g. Gene_Biotype = list( c(), c("coding", "processed_pseudogene" )).

This will allocate NULL/default settings to Intelligence_Sum_Stats and stratify the plot by only "coding" and "processed_pseudogene" regions for Household_Income_Sum_Stats.

Please note that HG19 is passed as the genome build for Intelligence_Sum_Stats but HG38 for Household_Income_Sum_Stats due to differences in the data sets' references.

```{r, results="hide", warning = FALSE, message = TRUE}

  setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

  Regional_Plots <- Regional_Plot(Data = c(Intelligence_Sum_Stats, Household_Income_Sum_Stats),  Chromosomes = c(9),  Genome_Build = c("grch37", "grch38"), 
                                  Recombination_Line_Colour = c("darkgrey", "pink"),  Point_Colour = c("blue", "red"), Y_Axis_Title = c("SigVal"),
                                  Gene_Biotype = list( c(), c("coding", "processed_pseudogene" )), Gene_Biotype_Colour = list( c(), c("blue", "green")) )
                                  

  setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

  Regional_Plots <- do.call(c, Regional_Plots)
  
  for (plot_name in names(Regional_Plots)) {
  plot <- Regional_Plots[[plot_name]]
  file_safe <- gsub("[^A-Za-z0-9]", "_", plot_name)
  filename <- paste0(file_safe, "_LmultiD", ".jpg")
  height_to_use <- attr(plot, "dynamic_height")
  if (is.null(height_to_use)) height_to_use <- 30

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = 30,
    height = height_to_use,   #nice divide by exon size to get ratio
    units = "in",
    dpi = 100,
    limitsize = FALSE
  )
  
}

```
Let's have a look at both chromosome 9 plots for our exemplar data sets.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

  setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")
  
  img_dir <- "C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots"
  
  # List files that end with 'Regional_Plot.jpg'
  img_paths <- list.files(
    path = img_dir,
    pattern = "Regional_Plot_LmultiD\\.jpg$",  # exact match
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Sort by modification time (oldest first)
  img_paths <- img_paths[order(file.info(img_paths)$mtime)]
  
  # Display all matching plots
  knitr::include_graphics(img_paths)
  
```
Single_Plot() also accepts a similar logic when creating multiple plots with separate customisation simultaneously.

```{r, results="hide", warning = FALSE, message = TRUE}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

Manhattan_Plot <- Single_Plot(Data = Sum_Stats, X_Axis_Title = "Region", Label_Height = 15, Anchor_Label = "centre", 
                              Label_Angle = c(90,45), Colour_Of_Index = "blue",  Chromosome_Colours = list( c("red", "orange"), c("purple", "orange")  ))

output_files <- c()
for (plot_name in names(Manhattan_Plot)) {
  file_path <- paste0(plot_name, "Complex_Single", ".jpg")
  ggplot2::ggsave(
    filename = file_path,
    plot = Manhattan_Plot[[plot_name]],
    width = 30, height = 15, units = "in", dpi = 100
  )
  output_files <- c(output_files, file_path)
}

```
```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")
knitr::include_graphics(output_files)

```

We can also use the MiamiR package to create a Miami plot (mirrored, paired Manhattan Plot) using two sets of summary statistics using the Miami_Plot() function.

We will be using the same pair of data as above; Intelligence_Sum_Stats and Household_Income_Sum_Stats, included in the MiamiR package.

This function automatically aligns and buffers coordinates across different data sets and genomic builds to create a neat visualisation. Let's look at a default plot, specifying both data sets where needed.

```{r, results="hide", warning = FALSE, message = TRUE}

Miami_Plot <- Miami_Plot(Top_Data = Intelligence_Sum_Stats, Bottom_Data = Household_Income_Sum_Stats)

setwd( "C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

ggplot2::ggsave("Test_Miami.jpg", plot = Miami_Plot, width = 30,
                height = 20, units = "in", dpi = 100)
                
```
Let's inspect our Miami Plot.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")
knitr::include_graphics("Test_Miami.jpg")

```
Miami_Plot(), similarly to Regional_Plot() is based on and inherits the functionality of Single_Plot() and it can be customised using the exact same applicable arguments as previously seen.

However, it also allows for segregated and overall inheritance. If an inherited argument like Point_Size is modified, then it will apply to the top and bottom plots within the Miami Plot.

If the argument is prefixed with Top_ or Bottom_ then these customisations will be allocated accordingly. Let's start with overall customisation.

```{r, results="hide", warning = FALSE, message = TRUE}

Miami_Plot <- Miami_Plot(Top_Data = Intelligence_Sum_Stats, Bottom_Data = Household_Income_Sum_Stats, Chromosome_Colours = c("limegreen", "darkorange"),
                         Chromosome_Label_Drops = c(15,19,22), Diamond_Index = TRUE, Chromosome_Label_Size = 12)

setwd( "C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

ggplot2::ggsave("Custom_Test_Miami.jpg", plot = Miami_Plot, width = 30,
                height = 20, units = "in", dpi = 100)
                
```
Let's inspect our customised Miami Plot.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")
knitr::include_graphics("Custom_Test_Miami.jpg")

```
Now let's use built in argument prefixing to modify the top and bottom plots independently within the Miami_Plot() function, whilst also applying global customisation.

```{r, results="hide", warning = FALSE, message = TRUE}

Miami_Plot <- Miami_Plot(Top_Data = Sum_Stats[1], Bottom_Data = Sum_Stats[2], Top_Chromosome_Colours = c("purple", "grey"),
                         Top_Diamond_Index = TRUE, Bottom_Chromosome_Colours = c("pink", "red"), Chromosome_Label_Size = 30, Label_Size = 7 )

setwd( "C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

ggplot2::ggsave("Split_Custom_Test_Miami.jpg", plot = Miami_Plot, width = 30,
                height = 20, units = "in", dpi = 100)  
```
Let's inspect our split customisation within our new Miami Plot.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")
knitr::include_graphics("Split_Custom_Test_Miami.jpg")

```
The MiamiR package can also be used to inspect key SNPs in single or multiple sets of GWAS summary statistics by using the Forest_Plot() function. This function carefully aligns and manages SNP differences automatically and allows for much customisation.

Here's how to create a Forest Plot centered around effect values of key SNPs on chromosome 4 from the previously described Intelligence_Sum_Stats data set.

The function will decipher the test statistic type and defaults to ascertainment of the lead SNP in each independent region, which can be defined similarly to the arguments used in Regional_Plot().

```{r, results="hide", warning = FALSE, message = TRUE}

Forest_Plot <- Forest_Plot(Data = Intelligence_Sum_Stats, Chromosomes = 4)

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")   

ggplot2::ggsave(
  filename = "Initial_Forest.jpg",
  plot = Forest_Plot,
  height = attr(Forest_Plot, "dynamic_height", exact = TRUE),
  width = 10,
  dpi = 100,
  units = "in"
)

```
Let's inspect our Chromosome 4 Forest Plot.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")  

knitr::include_graphics("Initial_Forest.jpg")

```
We can also specify particular SNPs of interest to examine, in the order to be plotted. Let's look at specific, select SNPs across both Intelligence_Sum_Stats and Household_Income_Sum_Stats, whilst also demonstrating some of the customisation offered through MiamiR. Here we will use file paths, instead of dataframe objects to demonstrate the versatility of Forest_Plot().

```{r, results="hide", warning = FALSE, message = TRUE}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

Forest_Plot <- Forest_Plot(Data = Sum_Stats, order = c("P_Label", "BETA_CI_Label"), Right_Size = 12,
                           Selected_SNPs = c("rs2014952", "rs759550", "rs72832626", "rs6682851",                                              "rs142764139"), Data_Set_Colours = c("blue", "green"), Shapes = c("triangle", "circle"), 
                           Names = c("Intel", "Income"), Double_Label = T, Left_Fill = "lightblue")
                              
setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes")   

ggplot2::ggsave(
  filename = "Selected_Custom_Forest.jpg",
  plot = Forest_Plot,
  height = attr(Forest_Plot, "dynamic_height", exact = TRUE),
  width = 10,
  dpi = 100,
  units = "in"
)

```
Let's inspect our selective Forest_Plot.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")  

knitr::include_graphics("Selected_Custom_Forest.jpg")

```
The MiamiR package also allows for the same figures to be produced but for covariate effects from the raw outputs of base R statistical models, particularly useful for GWAS sensitivity analysis and epidemiological profiling. 

By specifying Model_Reference = T within Forest_Plot() and passing such model objects as the Data argument a helper function called Model_Munge() will be called internally and automatically process such data. 

Let's create two models using a dataset included in the MiamiR package called Fake_Demo_Data and pass them to the Forest_Plot() function whilst also demonstrating additional customisation.

```{r, results="hide", warning = FALSE, message = TRUE}

Model_One <- lm(BMI ~ Age + Sex + Ethnicity + Location, data = Fake_Demo_Data)
Model_Two <- lm(Dementia ~ Age + Sex + Ethnicity + Location, data = Fake_Demo_Data)

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

Forest_Plot <- Forest_Plot(Data = c(Model_One, Model_Two), order = c("P_Label", "CI_Label"),                                       Right_Size = 12, Left_Size = 13, Model_Reference = T, Axis_Buffer = 0.05,
               Names = c("BMI", "Dementia"), Right_Fill = "lightgreen", Left_Fill = "turquoise" ,
               Strips = F, Styles = c("bold", "underline"), Null_Line_Type = "solid")
                              
setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")   

ggplot2::ggsave(
  filename = "Model_Custom_Forest.jpg",
  plot = Forest_Plot,
  height = attr(Forest_Plot, "dynamic_height", exact = TRUE),
  width = 10,
  dpi = 100,
  units = "in"
)

```

Let's inspect our model based Forest_Plot.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

setwd("C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")  

knitr::include_graphics("Model_Custom_Forest.jpg")

```
Model_Munge() can also be used independently of Forest_Plot() on such Model Objects. Let's investigate its use with "Model_One".

```{r, results="hide", warning = FALSE, message = TRUE}

Model_Munged <- Model_Munge(Model_Object = "Model_One")

```
Let's inspect our output from Model_Munge()

```{r, results="hide", warning = FALSE, message = FALSE}

print(head(Model_Munged))

```

Another helper function which is integrated in to the Single_Plot() function under the Auto_Lab setting is Annotate_Data() which uses the coordinates of your data and enables annotation with RSIDs for the index SNPs on each chromosome - these will be added to the data frame in a new column called Lab. This is particularly useful when creating plots from only coordinate files.

This queries an API in a similar way to the Auto_LD functionality in Reigonal_Plot() and as such it is important to bear in mind time constraints and requisite reference builds.

Here is how to use Labelled_Data in isolation.

```{r, results="hide", warning = FALSE, message = TRUE}

Annotated_Data <- Annotate_Data(Data = Intelligence_Sum_Stats, Genome_Build = "grch37")

```
Let's inspect our annotated data frame

```{r, results="hide", warning = FALSE, message = FALSE}

print(Annotated_Data[!is.na(Annotated_Data$Lab), c("CHR", "POS", "A1", "A2", "Lab")])

```

Last, but not least within the MiamiR package is the Tube_Alloys() function. This function combines Single_Plot(), Regional_Plot() and Forest_Plot() functionality in a blunt fashion to output overall plots, regions of interest and SNP and GWAS modelling information. The user must simply pass GWAS Summary Statistic(s) via any of the aforementioned routes in addition to either combined Phenotype/Covariate files or separate Phenotype and Covariate files if accessory analysis is desired. 

Here is how to run Tube_Alloys(), focusing on chromosome 9 on Intelligence_Sum_Stats & Household_Income_Sum_Stats in addition to supplying Fake_PHENOS_Binary & Fake_COVARS datasets included in MiamiR. The function will automatically apply defaults and auto-adjustments whilst applying focused arguments to only relevant functions where supplied data allows.

A special save helper function called save_all_plots is also included to facilitate easy output.

```{r, results="hide", warning = FALSE, message = TRUE}

Plots <- Tube_Alloys(Data = c(Intelligence_Sum_Stats, Household_Income_Sum_Stats), Phenos = Fake_PHENOS_Binary, 
                             Covars = Fake_COVARS, Chromosomes = 1)

setwd( "C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

Tube_Alloys <- invisible(save_all_plots(Plots, base_dir = "Tube_Alloys_Plots", dpi = 100))
                
```

Let's inspect all outputs.

```{r, results="hide", warning = FALSE, message = FALSE,  fig.width=30, fig.height=15, dpi=100, out.width="100%"}

setwd( "C:/Users/callumon/Miami_Package_R/MiamiR/vignettes/Plots")

# Flatten to one character vector
all_files <- unlist(Tube_Alloys, use.names = FALSE)

# Embed all images in your RMarkdown/Notebook
knitr::include_graphics(all_files)

```
A utility function also included in the MiamiR package is the METASOFT_File_Gen() function which takes multiple data sets as inputs and munges them by default to the METASOFT meta-analysis format (wrapper TBC). Other formats may be specified and allele matching is automatically taken care of, similarly to the Forest_Plot() function. Here is how to use the simple function using our usual data sets. 

```{r, results="hide", warning = FALSE, message = TRUE}

METASOFT <- METASOFT_File_Gen(Data = c("Intelligence_Sum_Stats", "Household_Income_Sum_Stats"))

```
Let's preview the output 

```{r, results="hide", warning = FALSE, message = FALSE}

print(na.omit(METASOFT))

```
