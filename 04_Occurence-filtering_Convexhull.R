####################### Occurrence filtering & Convexhull
library(terra)
library(openxlsx)
library(sf)
library(Rarity) # plot des correlations multiples
library(virtualspecies)
library(sp)
library(geodata)
library(ggplot2)
library(reshape2) 
library(gridExtra)
library(biomod2)
library(maps)
library(tidyterra)
library(ggtext)

# Load Environmental Space
envir.space <- readRDS("data/raw/bioclim/baseline/Environmental_Space.rds")
# Load variable intervals
intervals <- readRDS("data/raw/bioclim/baseline/intervals.rds")
# Load environmental raster stack
Rastack <- rast("data/raw/bioclim/baseline/final_baseline.tif")

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$x

# Change baseline Spat raster as raster
baseline_raster <- as(Rastack, "Raster")

# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  Occu <- read.xlsx(paste0("data/raw/occurences/Occurences_", Sp, ".xlsx")) # Load occurrences for this species
  
  ############# 1. Occurrence filtering #############
  
  ### a. Remove occurrences outside land
  
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
  
  ### b. Remove duplicated occurrences in Environmental Space
  
  # Extract environmental values at occurrence locations
  var.occ <- terra::extract(Rastack, Occu_2[, c("x", "y")])
  
  # Map occurrences to environmental intervals based on initial variable ranges
  # In which grid cells do the occurrences fall based on the initial intervals provided?
  # Where do the occurrences position themselves within the defined environmental space?
  env.itvl <- sapply(colnames(var.occ[,-1]), function(x, combs, seqs.)
    cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
    combs = var.occ[,-1], seqs. = var.intervals)
  
  duplicated.cells.itv <- which(duplicated(env.itvl)) # Identify duplicated environmental conditions
  
  message("\n\n ---- ", Sp, " ----\n",
          "Total number of occurrences with environmental conditions in the geographical space: ", 
          nrow(env.itvl),
          "\nNumber of duplicated conditions: ", 
          length(duplicated.cells.itv),
          "\nNumber of unique occurrences (environmental space): ",
          nrow(env.itvl[-duplicated.cells.itv, ]),
          "\n\nTotal number of presences: ",
          length(which(Occu_2$Observed == 1)),
          "\nNumber of duplicated presences: ", 
          length(which(Occu_2$Observed[duplicated.cells.itv] == 1)),
          "\nNumber of unique presences: ",
          length(which(Occu_2$Observed[-duplicated.cells.itv] == 1)),
          "\n\nTotal number of absences: ",
          length(which(Occu_2$Observed == 0)),
          "\nNumber of duplicated absences: ", 
          length(which(Occu_2$Observed[duplicated.cells.itv] == 0)),
          "\nNumber of unique absences: ",
          length(which(Occu_2$Observed[-duplicated.cells.itv] == 0)))
  
  # Filter out duplicated occurrences if any
  if(length(duplicated.cells.itv)>0){
    Occu_3 <- Occu_2[-duplicated.cells.itv, ] ## Remove duplicated occurrences
    env.itvl_2 <- env.itvl[-duplicated.cells.itv, ] ## new environmental interval without duplicated ones
    var.occ_2 <- var.occ[-duplicated.cells.itv, ]
  }else{
    Occu_3 <- Occu_2
    env.itvl_2 <- env.itvl
    var.occ_2 <- var.occ
  } 
  
  ############# 2. Convex Hull  #############
  
  ### a. Convex Hull preparation
  # Replace interval values by their mid points.
  cur.sp.pixels <- sapply(colnames(env.itvl_2), function(x, int.to.replace, replacing.values){
    replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
  }, int.to.replace = env.itvl_2, replacing.values = envir.space$detailed.intervals) ### valeurs variables correspondant aux occurences uniques
  
  df_cur.sp.pixels <- as.data.frame (cur.sp.pixels)
  # Save occurrences to an Excel file
  
  xlsx::write.xlsx(df_cur.sp.pixels, paste0("data/Filtered_occurences/cur.sp.pixels_", Sp, ".xlsx"), row.names = F)
  
  # Remove occurrences outliers for the convex hull calculations (2,5% high and 2,5% low extreme)
  convex.hull.interval = 0.05 
  # Select for each variable 2.5% of the lower extreme and 2.5% of the higher extreme
  outs <- lapply(colnames(var.occ_2), function(x, df)
  {
    qt <- quantile(df[, x], probs = c(0 + convex.hull.interval / 2,
                                      1 - convex.hull.interval / 2))
    return(which(df[, x] <= qt[1] | df[, x] >= qt[2]))
  }, df = var.occ_2)
  
  outs <- unique(unlist(outs))
  # Remove occurrences outliers from midpoint dataframe
  cur.sp.pixels.filt <- cur.sp.pixels[-outs,]
  
  ### b. Convex Hull creation
  
  if(nrow(cur.sp.pixels.filt) >= 6){  # If there is enough points we use occurrence without outliers
    # else we use all occurrences
    cursp.convhull <- geometry::convhulln(cur.sp.pixels.filt,
                                          # options = "Qt"
    ) 
  }else{
    cursp.convhull <- convhulln(cur.sp.pixels[Occu_3$Observed == 1, ]) 
  }
  
  # Checking environmental conditions that are inside the convex hull vs. outside
  cursp.inhull <- geometry::inhulln(cursp.convhull, 
                                    as.matrix(envir.space$unique.conditions.in.env))
  
  saveRDS(cursp.inhull, file = paste0("data/ConvexHull/", Sp, "/cursp.inhull.rds"))
  
  
  ## ConvexHull data Visualisation
  # plot(Rastack[[1]])
  # convhull_coord_values <- envir.space$coords[cursp.inhull,]
  # convhull_coord_values <- as.data.frame(convhull_coord_values)
  # points(convhull_coord_values$x, convhull_coord_values$y, pch = 20, col = "blue", cex = 0.2) ## convex_hull
  
  presencepixels <- apply(envir.space$unique.conditions.in.env, 1,paste, collapse = " ")   %in% 
    ## For each row, the env values are concatenated into a single string with spaces between them.
    #check if elements of one vector are present in another. 
    apply(cur.sp.pixels.filt, 1, paste, collapse = " ")
  #It returns a logical vector (TRUE or FALSE), indicating whether each element in the first vector is found in the second vector
  
  saveRDS(presencepixels, file = paste0("data/ConvexHull/", Sp, "/presencepixels.rds"))
  
  # Final data frame of occurrences and their respective variable values
  Fin_occ_var <- cbind(Occu_3, var.occ_2)
  Fin_occ_var <- Fin_occ_var %>%
    dplyr::select(-ID)
  
  # Save occurrences to an Excel file
  xlsx::write.xlsx(Fin_occ_var, paste0("data/Filtered_occurences/Occ&Var_", Sp, ".xlsx"), row.names = F)
  
  ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
  
  ############# 3. Data visualization
  
  ## 3.1 Final occurrences
  ## 3.1.1 Map
  
  # png(paste0("output/Filt_occurrences_plot/Occurence_plot_", Sp, ".png"), width = 800, height = 600)  # Adjust dimensions as needed
  # 
  # plot(Rastack[[1]])
  # points(Occu_2[ , c("x", "y")], pch = 20, cex = 0.5, col = "black")
  # points(Occu_3[ , c("x", "y")], pch = 20, cex = 0.5, col = "#F44336")
  # # Add a title
  # title(main = paste0("Final occurrence of ", Sp))
  # # Close the PNG device and save the image
  # dev.off()
  # 
  ## 3.1.2 2D visualisation
  ## a. Plot with 2 environmental variables
  # 
  # # Plot raw data (var.occ) as red dots
  # plot(var.occ$CHELSA_bio7, var.occ$CHELSA_bio5,
  #      xlab = "CHELSA_bio7",  # Label for the x-axis
  #      ylab = "CHELSA_bio5",  # Label for the y-axis
  #      pch = 19,  # Points as solid dots
  #      col = "black",  # Color of the raw data points (missing data)
  #      main = paste0(Sp," Raw and Filtered Occurence Data (data deleted in Red)"),
  #      xlim = range(c(var.occ$CHELSA_bio7, var.occ$CHELSA_bio7)),
  #      ylim = range(c(var.occ$CHELSA_bio5, var.occ$CHELSA_bio5)))
  # 
  # # Add filtered data (var.occ_3) as black dots
  # points(var.occ_2$CHELSA_bio7, var.occ_2$CHELSA_bio5,
  #        pch = 19,  # Points as solid dots
  #        col = "#f44336")  # Color of the filtered data points
  # 
  # 
  # # # b. 2D Environmental space and convexhull
  # var.occ_convhull <-  var.occ_2[-outs,]
  # env_space <- data.frame(envir.space$unique.conditions.in.env)
  # 
  # ggplot() +
  #   #   # Add a rectangle for the extent of the env_space dataframe
  #   geom_rect(aes(xmin = min(env_space$CHELSA_bio5), xmax = max(env_space$CHELSA_bio5), 
  #                 ymin = min(env_space$CHELSA_bio7), ymax = max(env_space$CHELSA_bio7)),
  #             fill = "green", alpha = 0.5) +
  #   # Add a rectangle for the extent of the Convex_hull dataframe
  #   geom_rect(aes(xmin = min(var.occ_convhull$CHELSA_bio5), xmax = max(var.occ_convhull$CHELSA_bio5), 
  #                 ymin = min(var.occ_convhull$CHELSA_bio7), ymax = max(var.occ_convhull$CHELSA_bio7)),
  #             fill = "blue") +
  #   # Plot the points for occurences as black dots
  #   geom_point(data = var.occ_2, aes(x = CHELSA_bio5, y = CHELSA_bio7), color = "black", size = 0.5) +
  #   # Plot points for the Convex_hull dataframe
  #   geom_point(data = var.occ_convhull, aes(x = CHELSA_bio5, y = CHELSA_bio7), color = "red",size = 0.5) +
  #   #   # Customize plot limits for better visualization
  #   xlim(min(c(env_space$CHELSA_bio5, env_space$CHELSA_bio5)) - 1, 
  #        max(c(env_space$CHELSA_bio5, env_space$CHELSA_bio5)) + 1) +
  #   ylim(min(c(env_space$CHELSA_bio7, env_space$CHELSA_bio7)) - 1, 
  #        max(c(env_space$CHELSA_bio7, env_space$CHELSA_bio7)) + 1) +
  #   #   # Add labels and title
  #   labs(title = "Extent Comparison Between env_space and Convex_hull",
  #        x = "CHELSA_bio5", y = "CHELSA_bio7") +
  #   theme_minimal()
  # 
  ## 3.2 Environmental range
  
  # Reshape the data into long format for ggplot
  var.occ_2_long <- var.occ_2[, -1] %>%
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
  # 
  # Loop over each environmental variable and create a boxplot for each
  for (var in unique(var.occ_2_long$Variable)) {
    plot_obj <- ggplot(var.occ_2_long[var.occ_2_long$Variable == var, ], aes(x = Variable, y = Values)) +
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
  # 
  # Arrange all the plots in a grid
  all_plots <- grid.arrange(grobs = plot_list, ncol = 3,
                            top = Sp)
  # # Save the grid to a PNG file
  ggsave(paste0("output/Filt_occurrences_plot/Variable_response_", Sp, ".png"), plot = all_plots, width = 10, height = 8, dpi = 300)
  # 
  ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
}
