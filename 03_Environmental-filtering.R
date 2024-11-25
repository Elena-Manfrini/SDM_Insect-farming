library(terra)
# library(sf)
# library(Rarity) # plot des correlations multiples
# library(virtualspecies)
# library(sp)
# library(geodata)
library(openxlsx)
library(ggplot2)
# library(reshape2) 
library(gridExtra)
library(biomod2)
# library(maps)
# library(tidyterra)
# library(ggtext)

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
    require(rgl)
    rgl::plot3d(all.env.pixels[, vars.to.plot])
  }
  # Return a list containing detailed intervals, unique conditions, and coordinates
  return(list(detailed.intervals = possible.combs,
              unique.conditions.in.env = all.env.pixels,
              coords = all.xy))
}

# Function to calculate midpoints of intervals
midpoints <- function(x) {
  lower <- as.numeric(gsub(",.*", "", gsub("\\(|\\[|\\)|\\]", "", x)))
  upper <- as.numeric(gsub(".*,", "", gsub("\\(|\\[|\\)|\\]", "", x)))
  midpoint <- lower + (upper - lower) / 2 ### ????? can't we just do the mean ??????
  return(midpoint)
}

#########################################################################################


############# 1. Environmental Space
# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

# Get values from the raster stack, including NA for missing values (e.g., ocean areas)
values <- values(Rastack) 
combinations <- as.data.frame(values)
combinations <- combinations[complete.cases(combinations), ] # Remove rows with NA


step <- 0.02

# Calculate quantiles for each environmental variable to set up intervals
# Each variable's values are split into 100 equal parts (for 100 quantile intervals)
bio5 <- quantile(combinations[,1], probs=seq(0,1,step))
# plot(bio5)
bio7 <- quantile(combinations[,2], probs=seq(0,1,step))
hurs_min <- quantile(combinations[,3], probs=seq(0,1,step))
hurs_range <- quantile(combinations[,4], probs=seq(0,1,step))
npp <- quantile(combinations[,5], probs=seq(0,1,step))
croplands <- quantile(combinations[,6], probs=seq(0,1,step))
# plot(Npp)

# Create a list of unique intervals for each variable (6-dimensional environmental grid)
intervals <- list(
  bio5 = sort(unique(as.numeric(bio5))),
  bio7 = sort(unique(as.numeric(bio7))),
  hurs_min = sort(unique(as.numeric(hurs_min))),
  hurs_range = sort(unique(as.numeric(hurs_range))),
  npp = sort(unique(as.numeric(npp))),
  croplands = sort(unique(as.numeric(croplands)))
)


# Assign names to intervals based on raster stack variable names
names(intervals) <- names(Rastack)
saveRDS(intervals, file = "data/intervals.rds")

# Calculate environmental combinations to create an environmental space
envir.space <- computeEnvCombinations(
  env.stack = Rastack,
  var.intervals = intervals, 
  plot = T,
  vars.to.plot = 1:5) 

saveRDS(envir.space, file = "data/Environmental_Space.rds")

############# 2. Occurrences filtering
Rastack <- rast("data/final_baseline.tif")
envir.space <- readRDS("data/Environmental_Space.rds")
intervals <- readRDS("data/intervals.rds")

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$x

# Change baseline Spat raster as raster
# baseline_raster <- as(Rastack, "Raster")
baseline_raster <- Rastack

# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  Occu <- read.xlsx(paste0("data/raw/occurences/Occurences_", Sp, ".xlsx")) # Load occurrences for this species

  ### 2.1 Remove occurrences outside land
  
  # Proceed only if there are more than 10 occurrences per environmental variable
  # if (length(Occu$x) > 60){
    # Rasterize occurrences to align with baseline raster
    Occu_r <- rasterize(as.matrix(Occu[,1:2]), baseline_raster)
    Occu_r <- as.data.frame(Occu_r, xy=T, na.rm = FALSE) # Convert to data frame with coordinates

    # Combine occurrences with the baseline raster values
    Rastab <- as.data.frame(baseline_raster[[1]], xy=T, na.rm = FALSE) # Take one layer of baseline raster
    Occu_r <- cbind(Occu_r,Rastab[[3]]) # Add variable values
    Occu_r <- Occu_r[complete.cases(Occu_r), ] # Remove occurrences outside land
    Occu_r <- Occu_r[, -4]
    Occu_r[3] <- 1 # Assign presence = 1 for each occurrence
    colnames(Occu_r) <- c("x","y","Observed") # Rename columns
  # }
  
   # var.intervals <- intervals # Environmental intervals for categorizing occurrence data
  
  ### 2.2 Remove duplicated occurrences in Environmental Space
  
  # Extract environmental values at occurrence locations
  var.occ <- terra::extract(Rastack, Occu_r[, c("x", "y")])
  var.occ <- var.occ[complete.cases(var.occ), ]
  
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
    Occu_3 <- Occu_r[-duplicated.cells.itv, ] ## Remove duplicated occurrences
    env.itvl_2 <- env.itvl[-duplicated.cells.itv, ] ## new environmental interval without duplicated ones
    var.occ_2 <- var.occ[-duplicated.cells.itv, ]
  } else {
    Occu_3 <- Occu_r
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

  ## ConvexHull data Visualisation
  # plot(Rastack[[1]])
  # convhull_coord_values <- envir.space$coords[cursp.inhull,]
  # convhull_coord_values <- as.data.frame(convhull_coord_values)
  # points(convhull_coord_values$x, convhull_coord_values$y, pch = 20, col = "blue", cex = 0.2) ## convex_hull
  
  # La ligne suivante permet de savoir quels pixels de l'espace vnrionnemental
  # correspondent à des présences (hors outliers) de notre espèce
  presencepixels <- apply(envir.space$unique.conditions.in.env, 1, paste, collapse = " ")   %in% 
    ## For each row, the env values are concatenated into a single string with spaces between them.
    #check if elements of one vector are present in another. 
    apply(cur.sp.pixels, 1, paste, collapse = " ")
  #It returns a logical vector (TRUE or FALSE), indicating whether each element in the first vector is found in the second vector
  
  # Final dataframe of occurences and their respective variable values
    Fin_occ_var <- cbind(Occu_3, var.occ_2)
    Fin_occ_var <- Fin_occ_var %>%
      dplyr::select(-ID)

    # Save occurrences to an Excel file
    if(!dir.exists("data/Filtered_occurences")) {
      dir.create("data/Filtered_occurences,")
    }
    xlsx::write.xlsx(Fin_occ_var, paste0("data/Filtered_occurences/Occ&Var_", Sp, ".xlsx"), row.names = F)
    
    ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    
    ############# 4. Data visualization
    
    ## 4.1 Final occurrences
    ## 4.1.1 Map

    # png(paste0("output/Filt_occurrences_plot/Occurence_plot_", Sp, ".png"), width = 800, height = 600)  # Adjust dimensions as needed
    # 
    # plot(Rastack[[1]])
    # points(Occu_r[ , c("x", "y")], pch = 20, cex = 0.5, col = "black")
    # points(Occu_3[ , c("x", "y")], pch = 20, cex = 0.5, col = "#F44336")
    # # Add a title
    # title(main = paste0("Final occurrence of ", Sp))
    # # Close the PNG device and save the image
    # dev.off()
    # 
    ## 4.1.2 2D visualisation
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
    # # b. 2D Environmental space and convexhull
    var.occ_convhull <-  var.occ_2[-outs,]
    env_space <- data.frame(envir.space$unique.conditions.in.env)
    
    ggplot() +
    #   # Add a rectangle for the extent of the env_space dataframe
       geom_rect(aes(xmin = min(env_space$bio5), xmax = max(env_space$bio5), 
                     ymin = min(env_space$bio7), ymax = max(env_space$bio7)),
                 fill = "green", alpha = 0.5) +
       # Add a rectangle for the extent of the Convex_hull dataframe
       geom_rect(aes(xmin = min(var.occ_convhull$bio5), xmax = max(var.occ_convhull$bio5), 
                     ymin = min(var.occ_convhull$bio7), ymax = max(var.occ_convhull$bio7)),
                 fill = "blue") +
       # Plot the points for occurences as black dots
       geom_point(data = var.occ_2, aes(x = bio5, y = bio7), color = "black", size = 0.5) +
       # Plot points for the Convex_hull dataframe
       geom_point(data = var.occ_convhull, aes(x = bio5, y = bio7), color = "red",size = 0.5) +
    #   # Customize plot limits for better visualization
       xlim(min(c(env_space$bio5, env_space$bio5)) - 1, 
            max(c(env_space$bio5, env_space$bio5)) + 1) +
       ylim(min(c(env_space$bio7, env_space$bio7)) - 1, 
            max(c(env_space$bio7, env_space$bio7)) + 1) +
    #   # Add labels and title
       labs(title = "Extent Comparison Between env_space and Convex_hull",
            x = "bio5", y = "bio7") +
       theme_minimal()
    # 
    ## 4.2 Environmental range
    
    # Reshape the data into long format for ggplot
     var.occ_2_long <- var.occ_2[, -1] %>%
       tidyr::pivot_longer(cols = everything(), names_to = "Variable", values_to = "Values")
 
    # Create a list to store each ggplot
     plot_list <- list()
    # Define a list of colors for each variable
     colors <- c("bio5" = "brown1",   # Blue
                 "bio7" = "brown",   # Orange
                 "hurs_min" = "#2BDBCA", # Green
                 "hurs_range" = "#488680", # Red
                 "npp" = "green4",     # Purple
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

############# 5. Pseudo absence generation

    number.PA <- nrow(Fin_occ_var)
    runs.PA <-5
    
    cursp.rundata <- Fin_occ_var[,-c(1:2)] # Environmental conditions
    
    pseudoabs.biomod.table <- matrix(FALSE,
                                     nrow = nrow(Fin_occ_var) + runs.PA * number.PA,,
                                     ncol = runs.PA) # creation table vide qui contient 5 pre run de pseudo-absence : meme nmbre pseudo-absence que de presence
    colnames(pseudoabs.biomod.table) <- paste0("PA", 1:runs.PA)
    pseudoabs.biomod.table[1:nrow(Fin_occ_var), ] <- TRUE
    
    # Occurrences coordinates
    cursp.xy <- Occu_3[, c("x", "y")] 
    
    for(PA in 1:runs.PA){
      
      # Sampling pseudo-absences outside the convex hull
      cursp.pseudoabs <- sample(which(!cursp.inhull & !presencepixels), 
                                size = number.PA,
                                replace = FALSE) 
      
      cursp.rundata <- rbind.data.frame(
        cursp.rundata,
        data.frame(Observed = NA,
                   envir.space$unique.conditions.in.env[cursp.pseudoabs, ]))
      
      cursp.xy <- rbind(cursp.xy,
                        data.frame(xyFromCell(Rastack[[1]], 
                                              cursp.pseudoabs)))
      
      pseudoabs.biomod.table[(nrow(cur.sp.pixels) + 1 + (PA - 1) * number.PA):
                               (nrow(cur.sp.pixels) + PA * number.PA), PA] <- TRUE
      
      #  #Combine environmental of observed and pseudoabscence data
      # cursp.rundata <- rbind.data.frame(cursp.rundata,
      #                                   data.frame(Observed = NA,
      #                                              envir.space$unique.conditions.in.env[cursp.pseudoabs, ]))
      # 
      # pseudoabs.biomod.table[(nrow(cur.sp.pixels) + 1 + (PA - 1) * number.PA):
      #                          (nrow(cur.sp.pixels) + PA * number.PA), PA] <- TRUE
      
      # pseudoabs.biomod.table[, PA] <- TRUE ##
      
           
      # #Combine environmental of observed and pseudoabscence coordinates
      # pseudoabscoord <- data.frame(envir.space$coords[cursp.pseudoabs,])
      # curcoordsp <- rbind(coord,pseudoabscoord)
      # data.frame(xyFromCell(Rastack[[1]],cursp.pseudoabs))) ## Get coordinates of NA + thoose of observed values
      
      # duplicated <- which(duplicated(curcoordsp))
      # 
      # test <- rbind.data.frame(Fin_occ_var[,c(1:3)],
      #                          + data.frame(Observed = NA,
      #                                       +            envir.space$coords[cursp.pseudoabs, ]))
      # test_6579 <- test[6579,]
      # test_6579 <- test_6579[,c(1,2)]
      # 
      # # Identify rows that are duplicates of row 31
      # duplicate_rows <- apply(test[,c(1,2)], 1, function(row) all(row == test_6579, na.rm = TRUE))
      # # Get the indices of duplicate rows
      # duplicate_indices <- which(duplicate_rows)
      # 
      # test[c(4827,6579),]
      
      
      # # Plot the base raster map
      # plot(Rastack[[1]])
      # # # Plot points from 'test' on the map
      # points(convhull_coord_values$x, convhull_coord_values$y, pch = 20, col = "blue", cex = 0.2) ## convex_hull 
      # points(curcoordsp$x, curcoordsp$y, pch = 20, col = "black", cex = 0.7) ## Pseudoabsences + presence
      # points(coord$x, coord$y, pch = 20, col = "red", cex = 0.7) ## Occurences 
      # # Add a legend to distinguish the point sets
      # legend("bottomleft", legend = c("Convex hull", "Pseudoabsence", "Occurences"), col = c("blue", "black", "red"), pch = 10, pt.cex = 0.5)
    }
      
    coorxy <- cursp.xy
    occurrences  <- cursp.rundata
    PATable <- pseudoabs.biomod.table
    # curocc.obs <- cursp.xy[,"Observed"] # Observed vaues : 1 & NA
    # curenv <- cursp.xy[,-1] # Environmental variables
      
      # Chemin de sauvegarde
      save_dir <- paste0("models/", Sp)
      
      if(!dir.exists(save_dir)) {
        dir.create(save_dir)
      }
      
      # Formatage des données pour BIOMOD2
      run_data <- BIOMOD_FormatingData(
        resp.name = Sp, 
        resp.var = occurrences$Observed, 
        expl.var = occurrences[, -1],   # Variables prédictives (rasterstack propre à l'espèce)
        dir.name = save_dir,  # Dossier de stockage des modèles
        resp.xy = coorxy,     # Coordonnées xy des présences et pseudo-absences
        PA.strategy = 'user.defined',
        PA.user.table = PATable)
      
      saveRDS(run_data, file = paste0("models/", Sp, "/run_data.RDS"))
      
      
      ##############################################################################
      # BLOCK CROSS VALIDATION
      
      proj_names <- readRDS(paste0("models/", Sp, "/run_data.RDS"))
      #proj_names
      
      table_cv <- bm_CrossValidation(
        bm.format = proj_names,
        strategy = "kfold",
        k = 5,
        nb.rep = 1,
        balance = "both")
      
      calib_summary <-
        summary(proj_names, calib.lines =  table_cv) %>%
        filter(dataset == "calibration") ## separation jeu de données en 5
      
      proj_names@data.env.var
      proj_names@data.species
      proj_names@coord
      
      iwp <- (10^6)^(1 - occurrences$Observed) ## parametre biomod
      
      RF_param_list <- NULL
      XGBOOST_param_list <- NULL
      MAXNET_param_list <- NULL
      
      for (cvrun in 1:nrow(calib_summary)) {
        prNum <- calib_summary$Presences[cvrun]
        bgNum <- calib_summary$Pseudo_Absences[cvrun]
        wt <- ifelse(occurrences$Observed == 1, 1, prNum / bgNum) # weight = 1 if occurrences or prNum/bgNum if pseudoabsences
       
        # Ajustement des paramètres pour Random Forest
        if (prNum < bgNum) {
          RF_param_list[[paste0("_", calib_summary$PA[[cvrun]], "_", calib_summary$run[[cvrun]])]] <- list(
            ntree = 2000,  # Augmentation du nombre d'arbres
            mtry = floor(sqrt(ncol(Projection_Class))),  # Ajustement de mtry pour capturer plus de variables
            sampsize = c("0" = prNum, "1" = prNum),
            replace = TRUE
          )
        } else {
          RF_param_list[[paste0("_", calib_summary$PA[[cvrun]], "_", calib_summary$run[[cvrun]])]] <- list(
            ntree = 2000,  # Augmentation du nombre d'arbres
            mtry = floor(sqrt(ncol(Projection_Class))),  # Ajustement de mtry pour capturer plus de variables
            sampsize = c("0" = bgNum, "1" = bgNum),
            replace = TRUE
          )
        }
        # Ajustement des paramètres pour XGBOOST
        XGBOOST_param_list[[paste0("_", calib_summary$PA[[cvrun]], "_", calib_summary$run[[cvrun]])]] <- list(
          nrounds = 1000,  # Augmentation des itérations pour une meilleure convergence
          eta = 0.05,  # Augmentation du taux d'apprentissage pour accélérer la convergence
          max_depth = 7,  # Augmentation de la profondeur pour capturer plus d'interactions
          subsample = 0.9,  # Utilisation d'une plus grande fraction des données
          objective = "binary:logistic",
          gamma = 1,  # Introduction de la complexité pour éviter le surapprentissage
          colsample_bytree = 0.8,  # Augmentation de la fraction des caractéristiques à chaque arbre
          min_child_weight = 1,
          weight = wt,
          verbose = 0
        )
        
        # Ajustement des paramètres pour MAXNET
        MAXNET_param_list[[paste0("_", calib_summary$PA[[cvrun]], "_", calib_summary$run[[cvrun]])]] <- list(
          l1_regularizer = 0.1,  # Ajout d'une légère régularisation L1
          l2_regularizer = 0.1,  # Ajout d'une légère régularisation L2
          use_sgd = TRUE,
          set_heldout = 0,
          verbose = TRUE
        )
      }
      
      Param_Models <- list(
        RF.binary.randomForest.randomForest = RF_param_list,
        XGBOOST.binary.xgboost.xgboost = XGBOOST_param_list,
        MAXNET.binary.maxnet.maxnet = MAXNET_param_list
      )
      
      
      
      model_parameters <- bm_ModelingOptions(
        data.type = "binary",
        models = c('RF', 'XGBOOST','MAXNET'),
        strategy = "user.defined",
        user.base = "default",
        user.val = list(
          RF.binary.randomForest.randomForest = RF_param_list,
          XGBOOST.binary.xgboost.xgboost = XGBOOST_param_list,
          MAXNET.binary.maxnet.maxnet = MAXNET_param_list),
        bm.format = proj_names,
        calib.lines = table_cv) 
      
      
      # Modélisation
      myBiomodModelOut <- BIOMOD_Modeling(
        run_data,
        modeling.id = "1",
        models = c('RF', 'XGBOOST','MAXNET'),
        OPT.strategy = "user.defined",
        OPT.user = model_parameters,
        CV.strategy = 'user.defined',
        CV.user.table = table_cv,
        CV.do.full.models = FALSE,
        var.import = 10,
        metric.eval = c('BOYCE'),
        do.progress = TRUE
      )
      
      # Sauvegarder le résultat de la modélisation
      saveRDS(myBiomodModelOut, file = paste0(save_dir, VecSp[i], "_model_output.rds"))
      
      # Obtenir les évaluations des modèles
      evals <- get_evaluations(myBiomodModelOut)
      
      # Filtrer les modèles avec un indice de Boyce > 0.7
      Choosen_Model <- evals %>% filter(validation > 0.7)
      
      # Récupération des projections correspondantes
      selected_models <- Choosen_Model$full.name # Nom des modèles sélectionnés
      
      
      resp <- bm_PlotResponseCurves(bm.out = myBiomodModelOut,
                                    models.chosen = selected_models,
                                    # show.variables = cur_vars,
                                    fixed.var = "mean",
                                    data_species = occurrences$Observed,
                                    do.plot = FALSE,
                                    do.progress = FALSE)$tab
      
      colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")
      
      # setwd(str_c("C:/Users/ElenaPC/Documents/Quentin_M2/models",Epoque[j],"/",Method[l],sep=""))
      
      p<-ggplot() +
        geom_line(data = resp[grep('RF', resp$Model), ],
                  aes(x = Var.value, y = Response, group = Model), col = "darkgreen",alpha = 0.2) +
        geom_line(data = resp[grep('XGBOOST', resp$Model), ],
                  aes(x = Var.value, y = Response, group = Model), col = 'darkblue', alpha = 0.2) +
        geom_line(data = resp[grep('MAXNET', resp$Model), ], 
                  aes(x = Var.value, y = Response, group = Model), col = 'darkred', alpha = 0.2) +
        facet_wrap(~ Variable, scales = "free_x") +
        theme_bw() +
        ylim(0, 1.1) +
        xlab("Variable value")
      
      ggsave(str_c("Relation_Environnement",VecSp[k],".jpeg",sep=""),
             p,
             # "jpeg",
             dpi = 2000,
             bg = NULL,
             width = 11,
             height = 8.5,
             units = "in"
      )
}