library(terra)
library(sf)
library(Rarity) # plot des correlations multiples
library(virtualspecies)
library(sp)
library(geodata)
library(openxlsx)
library(ggplot2)
library(reshape2) 
library(gridExtra) 

###################### Define computeEnvCombinations function ######################
computeEnvCombinations <- function(env.stack,
                                   var.intervals,
                                   plot = TRUE,
                                   vars.to.plot = 1:3)
{
  # Convert the raster stack values to a data frame
  combinations <- as.data.frame(values(Rastack[[c(names(Rastack))]]))
  # Get coordinates for all cells in the raster stack and remove NA cells
  all.xy <- xyFromCell(Rastack, 1:ncell(Rastack))
  all.xy <- all.xy[- which(is.na(combinations[, 1])), ]
  combinations <- combinations[- which(is.na(combinations[, 1])), ]
  
  # All possible combination in the environment
  comb.cat <- sapply(colnames(combinations), function(x, combs, seqs.)
    cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
    combs = combinations, seqs. = var.intervals)      
  
  message("Total number of cells with environmental conditions in the geographical space: ", nrow(combinations),
          "\nNumber of duplicated conditions: ", length(which(duplicated(comb.cat))),
          "\nNumber of unique cells (environmental space): ", nrow(comb.cat[-which(duplicated(comb.cat)), ]))
  
  # Calculate midpoint values for each environmental variable interval
  possible.combs <- lapply(var.intervals, function(x)
    data.frame(interval = cut(x, x, right = FALSE, include.lowest = TRUE),
               mids = midpoints(cut(x, x, right = FALSE, include.lowest = TRUE))))
  
  # Assign midpoint values for each environmental variable in the spatial grid
  all.env.pixels <- sapply(colnames(comb.cat), function(x, int.to.replace, replacing.values)
  {
    replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
  }, int.to.replace = comb.cat, replacing.values = possible.combs)
  
  # Remove duplicate cells in environmental space
  duplicated.cells <- which(duplicated(all.env.pixels))
  if(length(duplicated.cells)){
    all.env.pixels <- all.env.pixels[-duplicated.cells, ]
    all.xy <- all.xy[-duplicated.cells, ]
  }
  
  # Optional 3D plot of the selected environmental variables
  if(plot)
  {
    plot3d(all.env.pixels[, vars.to.plot])
  }
  # Return a list containing detailed intervals, unique conditions, and coordinates
  return(list(detailed.intervals = possible.combs,
              unique.conditions.in.env = all.env.pixels,
              coords = all.xy))
}

# Function to calculate midpoints of intervals
midpoints <- function(x, dp=2){
  lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
  upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
  return(round(lower+(upper-lower)/2, dp))
}

#########################################################################################


############# 1. Environmental Space
# Load environmental raster stack
Rastack <- rast("data/raw/bioclim/baseline/final_baseline.tif")

# Get values from the raster stack, including NA for missing values (e.g., ocean areas)
values <- values(Rastack) 
combinations <- as.data.frame(values)
combinations <- combinations[complete.cases(combinations), ] # Remove rows with NA

# Calculate quantiles for each environmental variable to set up intervals
# Each variable's values are split into 100 equal parts (for 100 quantile intervals)
Bio5 <- quantile(combinations[,1], probs=seq(0,1,0.01))
Bio7 <- quantile(combinations[,2], probs=seq(0,1,0.01))
Hurs_min <- quantile(combinations[,3], probs=seq(0,1,0.01))
Hurs_range <- quantile(combinations[,4], probs=seq(0,1,0.01))
Npp <- quantile(combinations[,5], probs=seq(0,1,0.01))
croplands <- quantile(combinations[,6], probs=seq(0,1,0.01))

# Create a list of unique intervals for each variable (5-dimensional environmental grid)
intervals <- list(
  Bio5 = sort(unique(as.numeric(Bio5))),
  Bio7 = sort(unique(as.numeric(Bio7))),
  Hurs_min = sort(unique(as.numeric(Hurs_min))),
  Hurs_range = sort(unique(as.numeric(Hurs_range))),
  Npp = sort(unique(as.numeric(Npp))),
  croplands = sort(unique(as.numeric(croplands)))
)

# Assign names to intervals based on raster stack variable names
names(intervals) <- names(Rastack)

# Calculate environmental combinations to create an environmental space
envir.space <- computeEnvCombinations(
  env.stack = Rastack[[c(names(Rastack))]],
  var.intervals = intervals, #grille definie
  plot = F,
  vars.to.plot = 1:5) 


############# 2. Occurrences data

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$x

# change baseline Spat raster as raster
baseline_raster <- as(Rastack, "Raster")

# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  Occu <- read.xlsx(paste0("data/raw/occurences/Occurences_", Sp, ".xlsx")) # Load occurrences for this species
  
  # Proceed only if there are more than 10 occurrences per environmental variable
  if (length(Occu$x) > 60){
    # Rasterize occurrences to align with baseline raster
    Occu <- rasterize(Occu[,1:2],baseline_raster)
    Occu <- as.data.frame(Occu, xy=T) # Convert to data frame with coordinates
    
    # Combine occurrences with the baseline raster values
    Rastab <- as.data.frame(baseline_raster[[1]],xy=T) # Take one layer of baseline raster
    Occu <- cbind(Occu,Rastab[[3]]) # Add variable values
    Occu <- Occu[complete.cases(Occu[ , c(colnames(Occu))]), ] # Remove occurrences outside land
    Occu <- Occu[,-c(4)] 
    Occu$layer <- rep(1,length(Occu$x)) # Assign presence = 1 for each occurrence
    colnames(Occu) <- c("x","y","Observed") # Rename columns
  }
  
  var.intervals <- intervals # Environmental intervals for categorizing occurrence data
  species.occurrence <- Occu
  
  # Extract environmental values at occurrence locations
  cur.occ <- extract(Rastack, species.occurrence[, c("x", "y")])
  
  # Map occurrences to environmental intervals based on initial variable ranges
  # In which grid cells do the occurrences fall based on the initial intervals provided?
  # Where do the occurrences position themselves within the defined environmental space?
  env.comb.cursp <- sapply(colnames(cur.occ), function(x, combs, seqs.)
    cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
    combs = cur.occ, seqs. = var.intervals)
  
  duplicated.cells <- which(duplicated(env.comb.cursp)) # Identify duplicated environmental conditions
  
  message("\n\n ---- ", Sp, " ----\n",
          "Total number of occurrences with environmental conditions in the geographical space: ", 
          nrow(env.comb.cursp),
          "\nNumber of duplicated conditions: ", 
          length(duplicated.cells),
          "\nNumber of unique occurrences (environmental space): ",
          nrow(env.comb.cursp[-duplicated.cells, ]),
          "\n\nTotal number of presences: ",
          length(which(species.occurrence$Observed == 1)),
          "\nNumber of duplicated presences: ", 
          length(which(species.occurrence$Observed[duplicated.cells] == 1)),
          "\nNumber of unique presences: ",
          length(which(species.occurrence$Observed[-duplicated.cells] == 1)),
          "\n\nTotal number of absences: ",
          length(which(species.occurrence$Observed == 0)),
          "\nNumber of duplicated absences: ", 
          length(which(species.occurrence$Observed[duplicated.cells] == 0)),
          "\nNumber of unique absences: ",
          length(which(species.occurrence$Observed[-duplicated.cells] == 0)))
  
  # Filter out duplicated occurrences if any
  if(length(duplicated.cells)>0){
    species.occ.filt <- species.occurrence[-duplicated.cells, ] ## Remove duplicated occurences
    env.comb.cursp <- env.comb.cursp[-duplicated.cells, ]
  }else{
    species.occ.filt<-species.occurrence 
  } 

  # Replace interval ranges by midpoints for unique occurrences
  cur.sp.pixels <- sapply(colnames(env.comb.cursp), function(x, int.to.replace, replacing.values){
    replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
  }, int.to.replace = env.comb.cursp, replacing.values = envir.space$detailed.intervals) ### valeurs variables correspondant aux occurences uniques
  
  cur.pres <- cur.occ[which(species.occurrence$Observed == 1), ]
 
   # Removing outliers for the convex hull calculations
  
  convex.hull.interval = 0.05 # pour eliminer les 2,5% des extremes hauts et bas
  
  outs <- lapply(colnames(cur.occ), function(x, df)
  {
    qt <- quantile(df[, x], probs = c(0 + convex.hull.interval / 2,
                                      1 - convex.hull.interval / 2))
    return(which(df[, x] <= qt[1] | df[, x] >= qt[2]))
  }, df = cur.pres)
  
  outs <- unique(unlist(outs)) ## cellules resultantes uniques les plus extremes
  cur.pres.filt <- cur.pres[-outs, ] ## retire les valeurs les plus extremes
  
  ##### 2.1 Visualize data ####
  # 2.1.a Datavisualisation of deleted occurences data
  # Plot raw data (cur.occ) as red dots
  
  plot(cur.occ$CHELSA_bio7, cur.occ$CHELSA_bio5,
       xlab = "CHELSA_bio7",  # Label for the x-axis
       ylab = "CHELSA_bio5",  # Label for the y-axis
       pch = 19,  # Points as solid dots
       col = "red",  # Color of the raw data points (missing data)
       main = paste0(Sp," Raw and Filtered Occurence Data (data deleted in Red)"),
       xlim = range(c(cur.occ$CHELSA_bio7, cur.pres.filt$CHELSA_bio7)),
       ylim = range(c(cur.occ$CHELSA_bio5, cur.pres.filt$CHELSA_bio5)))
  
  # Add filtered data (cur.pres.filt) as black dots
  points(cur.pres.filt$CHELSA_bio7, cur.pres.filt$CHELSA_bio5,
         pch = 19,  # Points as solid dots
         col = "black")  # Color of the filtered data points
  
  # 2.1.b Save a final dataset with coordinates and variable values
  species.occurrence$ID <- seq(1, nrow(species.occurrence))
  merged.occ <- merge(species.occurrence,cur.pres.filt , by = "ID")
  merged.occ <- merged.occ %>%
    select(-ID, -Observed)
  
  # Save occurrences to an Excel file
  xlsx::write.xlsx(merged.occ, paste0("data/Filtered_occurences/Occ&Var_", Sp, ".xlsx"), row.names = F)

  #  2.1.c Visualize environmental value range
  # Reshape the data into long format for ggplot
  env_vars_long <- melt(cur.pres.filt[,-1])
  
  # Create a list to store each ggplot
  plots <- list()
  # Define a list of colors for each variable
  colors <- c("CHELSA_bio5" = "brown1",   # Blue
              "CHELSA_bio7" = "brown",   # Orange
              "CHELSA_hurs_min" = "#2BDBCA", # Green
              "CHELSA_hurs_range" = "#488680", # Red
              "CHELSA_npp" = "green4",     # Purple
              "globalCropland_2010CE" = "#DB792C")  # Brown
  
  # Loop over each environmental variable and create a boxplot for each
  for (var in unique(env_vars_long$variable)) {
    p <- ggplot(env_vars_long[env_vars_long$variable == var, ], aes(x = variable, y = value)) +
      geom_boxplot(fill = colors[var], color = "black", outlier.shape = 16) +
      labs(x = NULL, y = NULL, title = var) +
      theme_minimal() +
      theme(axis.text.x = element_blank(),  # Remove x-axis labels for individual plots
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 10),  # Adjust size of y-axis labels
            axis.title.y = element_text(size = 12)) +
      scale_color_manual(values = colors)   +    # Use colors as specified
      guides(color = "none")  # Remove the color legend
    # Add the plot to the list
    plots[[var]] <- p
  }
  
  # Arrange all the plots in a grid
  all_plots <- grid.arrange(grobs = plots, ncol = 3,
                        top = Sp)
  # Save the grid to a PNG file
  ggsave(paste0("output/Filt_occurences_ggplot/Variable_response_", Sp, ".png"), plot = all_plots, width = 10, height = 8, dpi = 300)
  
  ##### 2.2 Remove duplicated occurence data in environmental variable intervals  ####
  env.comb.filt <- sapply(colnames(cur.pres.filt), function(x, combs, seqs.)
    cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
    combs = cur.pres.filt, seqs. = var.intervals) #### Dans quelles cases se situent les occurences (moins les extremes)
 
   # Remove duplicates
  duplicated.cells <- which(duplicated(env.comb.filt)) ### recreation des intervales = nouvelles cellules dupliquées.
  
  if(length(duplicated.cells)){
    env.comb.filt <- env.comb.filt[-duplicated.cells, ]
  } 
  
  # Replace values by mid points of intervals : chaque occurrences va avoir une valeur qui sera la distance par rapport a la mediane de l'hyper volume.
  cur.sp.pixels.filt <- sapply(colnames(env.comb.filt), function(x, int.to.replace, replacing.values)
  {
    replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
  }, int.to.replace = env.comb.filt, replacing.values = envir.space$detailed.intervals)
 
  
   # Calculating the convex hull of environmental conditions : considere une surface plus grande que les point d'occurence lorsque les pts d'occ sont trop rapprochés.
  if(nrow(env.comb.filt) >= 5){  # If there is enough points we use occurrence without outliers
    # else we use all occurrences
    cursp.convhull <- convhulln(cur.sp.pixels.filt,
                                options = "Qt"
    ) 
  }else{
    cursp.convhull <- convhulln(cur.sp.pixels[species.occ.filt$Observed == 1, ]) 
  }
  
  # Checking environmental conditions that are inside the convex hull vs. outside
  cursp.inhull <- inhulln(cursp.convhull, 
                          as.matrix(envir.space$unique.conditions.in.env))
  
  presencepixels <- apply(envir.space$unique.conditions.in.env, 1, 
                          paste, collapse = " ") %in%
    apply(cur.sp.pixels, 1, paste, collapse = " ")
}
