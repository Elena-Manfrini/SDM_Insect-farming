# ####################### Occurrence filtering & Convexhull
rm(list=ls())
library(terra)
library(openxlsx)
library(raster)
library(tidyr)
library(ggplot2)
library(gridExtra)

# Load Environmental Space
envir.space <- readRDS("data/Environmental_Space.rds")
# Load variable intervals
intervals <- readRDS("data/intervals.rds")
# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

# Change baseline Spat raster as raster
baseline_raster <- as(Rastack, "Raster")
Rastab <- as.data.frame(baseline_raster, xy=T, na.rm = FALSE) # Take one layer of baseline raster

i <- 1
# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  Occu <- read.xlsx(paste0("data/raw/occurences/Occurences_", Sp, ".xlsx")) # Load occurrences for this species
  ############# 2. Occurrences filtering
    ### 2.1 Remove occurrences outside land
    
    # Proceed only if there are more than 10 occurrences per environmental variable
    if (length(Occu$x) > 60){
    # Rasterize occurrences to align with baseline raster
    Occu_r<- rasterize(as.matrix(Occu[,1:2]), baseline_raster)
    Occu_r <- as.data.frame(Occu_r, xy=T, na.rm = FALSE) # Convert to data frame with coordinates
    
    # Combine occurrences with the baseline raster values
    # Rastab <- as.data.frame(baseline_raster[[4]], xy=T, na.rm = FALSE) # Take one layer of baseline raster
    Occu_r <- cbind(Occu_r,Rastab) # Add variable values
    Occu_r <- Occu_r[complete.cases(Occu_r), ] # Remove occurrences outside land
    Occu_r <- Occu_r[, c("x", "y", "layer")]
    Occu_r[3] <- 1 # Assign presence = 1 for each occurrence
    colnames(Occu_r) <- c("x","y","Observed") # Rename columns
    
    ### 2.2 Remove duplicated occurrences in Environmental Space
    
    # Extract environmental values at occurrence locations
    var.occ <- terra::extract(Rastack, Occu_r[, c("x", "y")])
    
    # Map occurrences to environmental intervals based on initial variable ranges
    # In which grid cells do the occurrences fall based on the initial intervals provided?
    # Where do the occurrences position themselves within the defined environmental space?
    env.itvl <- sapply(colnames(var.occ[,-1]), function(x, combs, seqs.)
      cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
      combs = var.occ[,-1], seqs. = intervals)
    
    duplicated.cells.itv <- which(duplicated(env.itvl)) # Identify duplicated environmental conditions
    
    message("\n\n ---- ", Sp, " ----\n",
            "Total number of occurrences with environmental conditions in the geographical space: ", 
            nrow(env.itvl),
            "\nNumber of duplicated conditions: ", 
            length(duplicated.cells.itv),
            "\nNumber of unique occurrences (environmental space): ",
            nrow(env.itvl[-duplicated.cells.itv, ]),
            "\n\nTotal number of presences (raw): ",
            nrow(Occu),
            "\n\nTotal number of presences (rasterized): ",
            length(which(Occu_r$Observed == 1)),
            "\nNumber of duplicated presences: ", 
            length(which(Occu_r$Observed[duplicated.cells.itv] == 1)),
            "\nNumber of unique presences: ",
            length(which(Occu_r$Observed[-duplicated.cells.itv] == 1)),
            "\n\nTotal number of absences: ",
            length(which(Occu_r$Observed == 0)),
            "\nNumber of duplicated absences: ", 
            length(which(Occu_r$Observed[duplicated.cells.itv] == 0)),
            "\nNumber of unique absences: ",
            length(which(Occu_r$Observed[-duplicated.cells.itv] == 0)))
    
    # Filter out duplicated occurrences if any
    if(length(duplicated.cells.itv) > 0){
      Occu_2 <- Occu_r[-duplicated.cells.itv, ] ## Remove duplicated occurrences
      env.itvl_2 <- env.itvl[-duplicated.cells.itv, ] ## new environmental interval without duplicated ones
      var.occ_2 <- var.occ[-duplicated.cells.itv, ]
    } else {
      Occu_2 <- Occu_r
      env.itvl_2 <- env.itvl
      var.occ_2 <- var.occ
    } 
    
    ############# 3. Convex Hull
    
    ### 3.1 Convex Hull preparation
    # Replace interval values by their mid points.
    cur.sp.pixels <- sapply(colnames(env.itvl_2), function(x, int.to.replace, replacing.values){
      replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
    }, int.to.replace = env.itvl_2, replacing.values = envir.space$detailed.intervals) ### valeurs variables correspondant aux occurences uniques
    
    # Remove occurrences outliers for the convex hull calculations (2,5% high and 2,5% low extreme)
    convex.hull.interval = 0.05
    
    Fin_occ_var <- cbind(Occu_2, var.occ_2)
    Fin_occ_var <- Fin_occ_var %>%
      dplyr::select(-ID)
    
    # Hurs
    hurs_min_q <- quantile(Fin_occ_var$hurs_min, probs = c(0.025, 0.975))
    hurs_min_ext <- Fin_occ_var[which(Fin_occ_var[, "hurs_min"] <= hurs_min_q[1] | Fin_occ_var[, "hurs_min"] >= hurs_min_q[2]),]
    
    # bio5
    bio5_q <- quantile(Fin_occ_var$bio5, probs = c(0.025, 0.975))
    bio5_ext <- Fin_occ_var[which(Fin_occ_var[, "bio5"] <= bio5_q[1] | Fin_occ_var[, "bio5"] >= bio5_q[2]),]
    
    # globalCropland_2010CE
    globalCropland_2010CE_q <- quantile(Fin_occ_var$globalCropland_2010CE, probs = c(0.025, 0.975))
    globalCropland_2010CE_ext <- Fin_occ_var[which(Fin_occ_var[, "globalCropland_2010CE"] <= globalCropland_2010CE_q[1] | Fin_occ_var[, "globalCropland_2010CE"] >= globalCropland_2010CE_q[2]),]
    
    # npp
    npp_q <- quantile(Fin_occ_var$npp, probs = c(0.025, 0.975))
    npp_ext <- Fin_occ_var[which(Fin_occ_var[, "npp"] <= npp_q[1] | Fin_occ_var[, "npp"] >= npp_q[2]),]
    
    # Plot
    
    # jpeg("output/InOut_convexhull_2,5per.jpeg", width = 15, height = 8.5, units = "in", res = 500)
    temp_raster <- Rastack[[1]]
    temp_raster[!is.na(temp_raster)] <- 1
    
    plot(temp_raster, col = "lightgrey", legend = FALSE)
    points(bio5_ext[, c("x", "y")], pch = 20, cex = 0.5, col = "red")
    points(hurs_min_ext[, c("x", "y")], pch = 20, cex = 0.5, col = "blue")
    points(globalCropland_2010CE_ext[, c("x", "y")], pch = 20, cex = 0.5, col = "orange")
    points(npp_ext[, c("x", "y")], pch = 20, cex = 0.5, col = "green")
    
    legend("left", legend = c("bio5", "hurs_min", "Cropland", "NPP"),
           col = c("red", "blue", "orange", "green"), pch = 20, cex = 0.8)
    
    
    # dev.off()
    
    
    
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
    
    
    ### 3.2 Convex Hull creation
    
    if(nrow(cur.sp.pixels.filt) >= 5){  # If there is enough points we use occurrence without outliers
      # else we use all occurrences
      cursp.convhull <- geometry::convhulln(cur.sp.pixels.filt,
                                            # options = "Qt"
      ) 
    }else{
      cursp.convhull <- convhulln(cur.sp.pixels[Occu_2$Observed == 1, ]) 
    }
    
    # Checking environmental conditions that are inside the convex hull vs. outside
    cursp.inhull <- geometry::inhulln(cursp.convhull, 
                                      as.matrix(envir.space$unique.conditions.in.env))
    
    
    # Save ConvexHull
    saveRDS(cursp.inhull, paste0("data/convexhull/", Sp, "_cursp.inhull.rds"))
    
    ## ConvexHull data Visualisation
    # plot(Rastack[[1]])
    # convhull_coord_values <- envir.space$coords[cursp.inhull,]
    # convhull_coord_values <- as.data.frame(convhull_coord_values)
    # points(convhull_coord_values$x, convhull_coord_values$y, pch = 20, col = "blue", cex = 0.2) ## convex_hull
    
    # Check where species occurences are located into the environmental space (using midpoints variable values)
    presencepixels <- apply(envir.space$unique.conditions.in.env, 1, paste, collapse = " ")   %in% 
      ## For each row, the env values are concatenated into a single string with spaces between them.
      # check if elements of one vector are present in another. 
      apply(cur.sp.pixels, 1, paste, collapse = " ")
    #It returns a logical vector (TRUE or FALSE), indicating whether each element in the first vector is found in the second vector
    
    # Save ConvexHull
    saveRDS(presencepixels, paste0("data/convexhull/", Sp, "_presencepixels.rds"))
    
    # Final dataframe of occurences and their respective variable values
    Fin_occ_var <- cbind(Occu_2, var.occ_2)
    Fin_occ_var <- Fin_occ_var %>%
      dplyr::select(-ID)
    
    outconv <- Fin_occ_var[outs,]
    nrow(outconv)
    inconv <- Fin_occ_var[-outs,]
    nrow(inconv)

    jpeg("output/InOut_convexhull_2,5per.jpeg", width = 15, height = 8.5, units = "in", res = 500)
    plot(Rastack[[1]])
    points(inconv[, c("x", "y")], pch = 20, cex = 0.5, col = "blue")
    points(outconv[, c("x", "y")], pch = 20, cex = 0.5, col = "red")
    dev.off()
    
    
    # Save occurrences to an Excel file
    if(!dir.exists("data/filtered_occurences")) {
      dir.create("data/filtered_occurences,")
    }
    xlsx::write.xlsx(Fin_occ_var, paste0("data/filtered_occurences/Occ&Var_final_15", Sp, ".xlsx"), row.names = F)
    
    ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    
    ############# 4. Data visualization
    
    ## 4.1 Final occurrences
    ## 4.1.1 Map
    
    # png(paste0("output/Filt_occurrences_plot/Occurence_plot_", Sp, ".png"), width = 800, height = 600)  # Adjust dimensions as needed
    # 
    # plot(Rastack[[1]])
    # points(Occu_r[ , c("x", "y")], pch = 20, cex = 0.5, col = "black")
    # points(Occu_2[ , c("x", "y")], pch = 20, cex = 0.5, col = "#F44336")
    # # Add a title
    # title(main = paste0("Final occurrence of ", Sp))
    # # Close the PNG device and save the image
    # dev.off()
    # 
    ## 4.1.2 2D visualisation
    ## a. Plot with 2 environmental variables
    # 
    # # Plot raw data (var.occ) as red dots
    # plot(var.occ$bio5, var.occ$hurs_min,
    #      xlab = "bio5",  # Label for the x-axis
    #      ylab = "hurs_min",  # Label for the y-axis
    #      pch = 19,  # Points as solid dots
    #      col = "black",  # Color of the raw data points (missing data)
    #      main = paste0(Sp," Raw and Filtered Occurence Data (data deleted in Red)"),
    #      xlim = range(c(var.occ$bio5, var.occ$bio5)),
    #      ylim = range(c(var.occ$hurs_min, var.occ$hurs_min)))
    # 
    # # Add filtered data (var.occ_3) as black dots
    # points(var.occ_2$bio5, var.occ_2$hurs_min,
    #        pch = 19,  # Points as solid dots
    #        col = "#f44336")  # Color of the filtered data points
    # 
    # 
    # # b. 2D Environmental space and convexhull
    # var.occ_convhull <-  var.occ_2[-outs,]
    # env_space <- data.frame(envir.space$unique.conditions.in.env)
    # 
    # ggplot() +
    #   #   # Add a rectangle for the extent of the env_space dataframe
    #   geom_rect(aes(xmin = min(env_space$bio5), xmax = max(env_space$bio5), 
    #                 ymin = min(env_space$hurs_min), ymax = max(env_space$hurs_min)),
    #             fill = "green", alpha = 0.5) +
    #   # Add a rectangle for the extent of the Convex_hull dataframe
    #   geom_rect(aes(xmin = min(var.occ_convhull$bio5), xmax = max(var.occ_convhull$bio5), 
    #                 ymin = min(var.occ_convhull$hurs_min), ymax = max(var.occ_convhull$hurs_min)),
    #             fill = "blue") +
    #   # Plot the points for occurences as black dots
    #   geom_point(data = var.occ_2, aes(x = bio5, y = hurs_min), color = "black", size = 0.5) +
    #   # Plot points for the Convex_hull dataframe
    #   geom_point(data = var.occ_convhull, aes(x = bio5, y = hurs_min), color = "red",size = 0.5) +
    #   #   # Customize plot limits for better visualization
    #   xlim(min(c(env_space$bio5, env_space$bio5)) - 1, 
    #        max(c(env_space$bio5, env_space$bio5)) + 1) +
    #   ylim(min(c(env_space$hurs_min, env_space$hurs_min)) - 1, 
    #        max(c(env_space$hurs_min, env_space$hurs_min)) + 1) +
    #   #   # Add labels and title
    #   labs(title = "Extent Comparison Between env_space and Convex_hull",
    #        x = "bio5", y = "bio7") +
    #   theme_minimal()
    # 
    ## 4.2 Environmental range
    
    # Reshape the data into long format for ggplot
    var.occ_2_long <- var.occ_2[, -1] %>%
      tidyr::pivot_longer(cols = everything(), names_to = "Variable", values_to = "Values")
    
    # Create a list to store each ggplot
    plot_list <- list()
    # Define a list of colors for each variable
    colors <- c("bio5" = "brown",   
                "hurs_min" = "#2BDBCA", 
                "npp" = "green4",
                "bio6" = "black",
                # "Human_pop_2000" = "#8494FF",
                # "Human_footprint" = "blue",
                "globalCropland_2010CE" = "#E68613")  
     
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
    
    } else {
      # Stop the script if the number of rows in Occu is not greater than 60
      stop("The number of occurrences is less than or equal to 60. Exiting the script.")
    }
}

