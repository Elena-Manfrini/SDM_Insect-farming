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
plot(Bio5)
Bio7 <- quantile(combinations[,2], probs=seq(0,1,0.01))
Hurs_min <- quantile(combinations[,3], probs=seq(0,1,0.01))
Hurs_range <- quantile(combinations[,4], probs=seq(0,1,0.01))
Npp <- quantile(combinations[,5], probs=seq(0,1,0.01))
croplands <- quantile(combinations[,6], probs=seq(0,1,0.01))
plot(Npp)

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

# Environmental space 2 dimensions 
env_space <- as.data.frame(envir.space$unique.conditions.in.env[,c("CHELSA_bio5", "CHELSA_bio7")])

############# 2. Occurrences data filtering

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$x

# change baseline Spat raster as raster
baseline_raster <- as(Rastack, "Raster")

# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  Occu <- read.xlsx(paste0("data/raw/occurences/Occurences_", Sp, ".xlsx")) # Load occurrences for this species
  
  ### 2.1 Remove occurences outside land
  
  # Proceed only if there are more than 10 occurrences per environmental variable
  if (length(Occu$x) > 60){
    # Rasterize occurrences to align with baseline raster
    Occu_2 <- rasterize(Occu[,1:2],baseline_raster)
    Occu_2 <- as.data.frame(Occu_2, xy=T) # Convert to data frame with coordinates
    
    # Combine occurrences with the baseline raster values
    Rastab <- as.data.frame(baseline_raster[[1]],xy=T) # Take one layer of baseline raster
    Occu_2 <- cbind(Occu_2,Rastab[[3]]) # Add variable values
    Occu_2 <- Occu_2[complete.cases(Occu_2[ , c(colnames(Occu_2))]), ] # Remove occurrences outside land
    Occu_2 <- Occu_2[,-c(4)] 
    Occu_2$layer <- rep(1,length(Occu_2$x)) # Assign presence = 1 for each occurrence
    colnames(Occu_2) <- c("x","y","Observed") # Rename columns
  }
  
  var.intervals <- intervals # Environmental intervals for categorizing occurrence data
  
  ### 2.2 Remove duplicated occurrences in environmental conditions
  
  # Extract environmental values at occurrence locations
  var.occ <- terra::extract(Rastack, Occu_2[, c("x", "y")])
  
  # Map occurrences to environmental intervals based on initial variable ranges
  # In which grid cells do the occurrences fall based on the initial intervals provided?
  # Where do the occurrences position themselves within the defined environmental space?
  env.itvl <- sapply(colnames(var.occ[,-1]), function(x, combs, seqs.)
    cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
    combs = var.occ[,-1], seqs. = var.intervals)
  
  duplicated.cells <- which(duplicated(env.itvl)) # Identify duplicated environmental conditions
  
  message("\n\n ---- ", Sp, " ----\n",
          "Total number of occurrences with environmental conditions in the geographical space: ", 
          nrow(env.itvl),
          "\nNumber of duplicated conditions: ", 
          length(duplicated.cells),
          "\nNumber of unique occurrences (environmental space): ",
          nrow(env.itvl[-duplicated.cells, ]),
          "\n\nTotal number of presences: ",
          length(which(Occu_2$Observed == 1)),
          "\nNumber of duplicated presences: ", 
          length(which(Occu_2$Observed[duplicated.cells] == 1)),
          "\nNumber of unique presences: ",
          length(which(Occu_2$Observed[-duplicated.cells] == 1)),
          "\n\nTotal number of absences: ",
          length(which(Occu_2$Observed == 0)),
          "\nNumber of duplicated absences: ", 
          length(which(Occu_2$Observed[duplicated.cells] == 0)),
          "\nNumber of unique absences: ",
          length(which(Occu_2$Observed[-duplicated.cells] == 0)))
  
  # Filter out duplicated occurrences if any
  if(length(duplicated.cells)>0){
    Occu_3 <- Occu_2[-duplicated.cells, ] ## Remove duplicated occurrences
    env.itvl_2 <- env.itvl[-duplicated.cells, ] ## new environmental interval without duplicated ones
  }else{
    Occu_3 <- Occu_2 
  } 

  ### 2.3 Removing occurrences outliers for the convex hull calculations (2,5% high and 2,5% low extreme)
  
  # Replace interval ranges by midpoints for unique occurrences in the environmental space
  cur.sp.pixels <- sapply(colnames(env.itvl_2), function(x, int.to.replace, replacing.values){
    replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
  }, int.to.replace = env.itvl_2, replacing.values = envir.space$detailed.intervals) ### valeurs variables correspondant aux occurences uniques
  
  var.occ_2 <- var.occ[which(Occu_3$Observed == 1), ]
 
   # Removing outliers for the convex hull calculations
  convex.hull.interval = 0.05 # pour eliminer les 2,5% des extremes hauts et bas
  
  outs <- lapply(colnames(var.occ_2), function(x, df)
  {
    qt <- quantile(df[, x], probs = c(0 + convex.hull.interval / 2,
                                      1 - convex.hull.interval / 2))
    return(which(df[, x] <= qt[1] | df[, x] >= qt[2]))
  }, df = var.occ_2)
  
  outs <- unique(unlist(outs)) ## cellules resultantes uniques les plus extremes
  
  var.occ_3 <- var.occ_2[-outs, ] ## Environmental variables withouth extreme
  Occu_4 <- Occu_3[-outs, ] ## Occurences withouth extreme
  
  ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
  # 2.1.1 Visualize data ####
  # 2.1.1.a Datavisualisation of deleted occurences data
  # Plot raw data (var.occ) as red dots
  plot(var.occ$CHELSA_bio7, var.occ$CHELSA_bio5,
       xlab = "CHELSA_bio7",  # Label for the x-axis
       ylab = "CHELSA_bio5",  # Label for the y-axis
       pch = 19,  # Points as solid dots
       col = "red",  # Color of the raw data points (missing data)
       main = paste0(Sp," Raw and Filtered Occurence Data (data deleted in Red)"),
       xlim = range(c(var.occ$CHELSA_bio7, var.occ_3$CHELSA_bio7)),
       ylim = range(c(var.occ$CHELSA_bio5, var.occ_3$CHELSA_bio5)))
  
  # Add filtered data (var.occ_3) as black dots
  points(var.occ_3$CHELSA_bio7, var.occ_3$CHELSA_bio5,
         pch = 19,  # Points as solid dots
         col = "black")  # Color of the filtered data points
  
  # 2.1.1.b Data visualisation of restricted environmental space for 2 variables
  ggplot() +
    # Add a rectangle for the extent of the env_space dataframe
    geom_rect(aes(xmin = min(env_space$CHELSA_bio5), xmax = max(env_space$CHELSA_bio5), 
                  ymin = min(env_space$CHELSA_bio7), ymax = max(env_space$CHELSA_bio7)),
              fill = "green", alpha = 0.5) +
    # Add a rectangle for the extent of the Convex_hull dataframe
    geom_rect(aes(xmin = min(var.occ_3$CHELSA_bio5), xmax = max(var.occ_3$CHELSA_bio5), 
                  ymin = min(var.occ_3$CHELSA_bio7), ymax = max(var.occ_3$CHELSA_bio7)),
              fill = "lightblue") +
    # Plot the points for var.occ as red dots
    geom_point(data = var.occ, aes(x = CHELSA_bio5, y = CHELSA_bio7), color = "red", size = 1) +
    # Plot points for the Convex_hull dataframe
    geom_point(data = var.occ_3, aes(x = CHELSA_bio5, y = CHELSA_bio7), color = "black") +
    
    # Customize plot limits for better visualization
    xlim(min(c(env_space$CHELSA_bio5, env_space$CHELSA_bio5)) - 1, 
         max(c(env_space$CHELSA_bio5, env_space$CHELSA_bio5)) + 1) +
    ylim(min(c(env_space$CHELSA_bio7, env_space$CHELSA_bio7)) - 1, 
         max(c(env_space$CHELSA_bio7, env_space$CHELSA_bio7)) + 1) +
    # Add labels and title
    labs(title = "Extent Comparison Between env_space and Convex_hull",
         x = "CHELSA_bio5", y = "CHELSA_bio7") +
    theme_minimal()
  
  # # 2.1.1.c Save a final dataset with coordinates and variable values
  # Occu_2$ID <- seq(1, nrow(Occu_2))
  # merged.occ <- merge(Occu_2,var.occ_3 , by = "ID")
  # merged.occ <- merged.occ %>%
  #   dplyr::select(-ID, -Observed)
  # 
  # # Save occurrences to an Excel file
  # xlsx::write.xlsx(merged.occ, paste0("data/Filtered_occurences/Occ&Var_", Sp, ".xlsx"), row.names = F)

  #  2.1.1.c Visualize environmental value range
  # Reshape the data into long format for ggplot
  var.occ_3_long <- var.occ_3[, -1] %>%
    tidyr::pivot_longer(cols = everything(), names_to = "Variable", values_to = "Values")
  
  # Create a list to store each ggplot
  plot_list <- list()
  # Define a list of colors for each variable
  colors <- c("CHELSA_bio5" = "brown1",   # Blue
              "CHELSA_bio7" = "brown",   # Orange
              "CHELSA_hurs_min" = "#2BDBCA", # Green
              "CHELSA_hurs_range" = "#488680", # Red
              "CHELSA_npp" = "green4",     # Purple
              "globalCropland_2010CE" = "#DB792C")  # Brown
  
  # Loop over each environmental variable and create a boxplot for each
  for (var in unique(var.occ_3_long$Variable)) {
    plot_obj <- ggplot(var.occ_3_long[var.occ_3_long$Variable == var, ], aes(x = Variable, y = Values)) +
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
    plot_list[[var]] <- plot_obj
  }
  
  # Arrange all the plots in a grid
  all_plots <- grid.arrange(grobs = plot_list, ncol = 3,
                        top = Sp)
  # Save the grid to a PNG file
  ggsave(paste0("output/Filt_occurences_ggplot/Variable_response_", Sp, ".png"), plot = all_plots, width = 10, height = 8, dpi = 300)
  
  ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
  
  ### 2.4 Remove duplicated occurrence data in environmental variable intervals
  # New calcul of intervals
  env.itvl_2 <- sapply(colnames(var.occ_3[,-1]), function(x, combs, seqs.)
    cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
    combs = var.occ_3, seqs. = var.intervals) #### Dans quelles cases se situent les occurences (moins les extremes)
  
  # Remove duplicates
  duplicated.cells <- which(duplicated(env.itvl_2)) ### recreate intervals = new duplicated cells
  
  if(length(duplicated.cells)){
    env.itvl_3 <- env.itvl_2[-duplicated.cells, ]
    var.occ_4 <- var.occ_3[-duplicated.cells, ]
    Occu_5 <- Occu_4[-duplicated.cells, ]
  } 
  
  # Replace values by mid points of intervals : chaque occurrences va avoir une valeur qui sera la distance par rapport a la mediane de l'hyper volume.
  cur.sp.pixels.filt <- sapply(colnames(env.itvl_3), function(x, int.to.replace, replacing.values)
  {
    replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
  }, int.to.replace = env.itvl_3, replacing.values = envir.space$detailed.intervals)
  
  
  if(nrow(cur.sp.pixels.filt) >= 6){  # If there is enough points we use occurrence without outliers
    # else we use all occurrences
    cursp.convhull <- geometry::convhulln(cur.sp.pixels.filt,
                                options = "Qt"
    ) 
  }else{
    cursp.convhull <- convhulln(cur.sp.pixels[Occu_3$Observed == 1, ]) 
  }
  
  # Checking environmental conditions that are inside the convex hull vs. outside
  cursp.inhull <- geometry::inhulln(cursp.convhull, 
                          as.matrix(envir.space$unique.conditions.in.env))
  
  presencepixels <- apply(envir.space$unique.conditions.in.env, 1,paste, collapse = " ")   %in% ## For each row, the env values are concatenated into a single string with spaces between them.
    #check if elements of one vector are present in another. 
    apply(cur.sp.pixels.filt, 1, paste, collapse = " ") # same here
  #It returns a logical vector (TRUE or FALSE), indicating whether each element in the first vector is found in the second vector
  # e.g 280 env condition at species occurrences are found inside the convex hull (all occurrences in the case of "Alphitobius diaperinus")

  
  # Final dataframe of occurences and their respective variable values
    Fin_occ_var <- cbind(Occu_5, var.occ_4)
    Fin_occ_var <- Fin_occ_var %>%
    dplyr::select(-ID)

    # Save occurrences to an Excel file
    xlsx::write.xlsx(Fin_occ_var, paste0("data/Filtered_occurences/Occ&Var_", Sp, ".xlsx"), row.names = F)
    
    
    # Save the map as a PNG
    # Save the grid to a PNG file
    png(paste0("output/Filt_occurences_ggplot/Occurence_plot_", Sp, ".png"), width = 800, height = 600)  # Adjust dimensions as needed
    
    map("world", xlim = range(Rastab$x), ylim = range(Rastab$y)) # Map occurences
    points(Occu[ , c("x", "y")], pch = 20, cex = 0.5, col = "red")
    points(Occu_5[ , c("x", "y")], pch = 20, cex = 0.5, col = "darkgreen")
    
    # Add a title
    title(main = paste0("Occurrence of", Sp))
    # Add a legend
    legend("left",                        # Position of the legend
           legend = c("Final", "Deleted"),  # Labels
           col = c("darkgreen", "red"),       # Colors matching the points
           pch = 20,                          # Point symbol used in the plot
           pt.cex = 1.5,                      # Point size in the legend
           cex = 1,                           # Text size in the legend
           bty = "n")  
    
    # Close the PNG device and save the image
    dev.off()
}
  #################### 3. Pseudo absence generation
    
    number.PA <- nrow(Fin_occ_var)
    runs.PA <-5
    
    pseudoabs.biomod.table <- matrix(FALSE,
                                   nrow = nrow(Fin_occ_var) + runs.PA * number.PA,
                                   ncol = runs.PA) # creation table vide qui contient 5 pre run de pseudo-absence : meme nmbre pseudo-absence que de presence
  colnames(pseudoabs.biomod.table) <- paste0("PA", 1:runs.PA)
  pseudoabs.biomod.table[1:nrow(Fin_occ_var), ] <- TRUE ## = presence
  
  # presence + Envirnomental variables
  rundata <- data.frame(
    Observed = 1, # Non-duplicated species occurrences
    var.occ_4[,-1])
  
  # Occurences coordinates
  coord <- Occu_5[, c("x", "y")] 
  
  plot = F
  for(PA in 1:runs.PA){
    try({
      # Sampling pseudo-absences outside the convex hull
      # AND not on presence points
      pseudoabs <- sample(which(!cursp.inhull & !presencepixels), 
                                size = number.PA,
                                replace = FALSE) #
      
      rundata <- rbind.data.frame(rundata,
                                        data.frame(Observed = NA,
                                                   envir.space$unique.conditions.in.env[pseudoabs, ])) # Bind environmental data of observed and randomly selected.
      
      coord <- rbind(coord,
                        data.frame(xyFromCell(Rastack[[1]],pseudoabs))) ## Get coordinates of NA + thoose of observed values
      
      pseudoabs.biomod.table[(nrow(cur.sp.pixels.filt) + 1 + (PA - 1) * number.PA):
                               (nrow(cur.sp.pixels.filt) + PA * number.PA), PA] <- TRUE
    })
    
     if(plot){
       plotconvexhull(allenvpixels = as.data.frame(envir.space$unique.conditions.in.env),
                      cursppixels = cur.sp.pixels.filt,
                      pseudoabs = pseudoabs,
                      curspinhull = cursp.inhull)
     }
    
    results <- list(occurrence.environment.matrix = rundata,
                    pseudoabs.biomod.table = pseudoabs.biomod.table, # pour chaque set pseudo abs : occurences + le set de pseudo abs
                    xy.coordinates = coord) ## creation liste avec coordonnées des occurences + pseudo abs (5 car 5 echantillonages) 
    
    coorxy <- results$xy.coordinates
    
    occurrences <- results$occurrence.environment.matrix
    
    PATable <- results$pseudoabs.biomod.table
    
  }
}
