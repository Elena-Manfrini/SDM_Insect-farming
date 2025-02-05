rm(list=ls())
library(terra)
library(openxlsx)
library(biomod2)
library(dplyr)

# Load Environmental Space
envir.space <- readRDS("data/Environmental_Space.rds")

# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

Rastab <- as.data.frame(Rastack, xy=T, na.rm = FALSE) # Take one layer of baseline raster

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

i <- 1
# Loop over each species to process occurrence data
 for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  
  # Occurrence and environmental values for the species
  Fin_occ_var <- read.xlsx(paste0("data/filtered_occurences/Occ&Var_final", Sp, ".xlsx"))
  
  # Convex hull and presence pixels data
  cursp.inhull <- readRDS(paste0("data/convexhull/", Sp, "_cursp.inhull.rds"))
  presencepixels <- readRDS(paste0("data/convexhull/", Sp, "_presencepixels.rds"))
  
  # Set the number of presence-absence (PA) points and the number of random runs for PA sampling
  number.PA <- nrow(Fin_occ_var) # Number of presence points = Number of pseudo absence points
  runs.PA <-5 # Number of random runs for pseudo-absence sampling
  
  # Prepare environmental data for the species 
  cursp.rundata <- Fin_occ_var[,-c(1:2)] # Environmental conditions
  
  # XY coordinates of the species' presence points
  cursp.xy <- Fin_occ_var[, c("x", "y")] 
  
  # Table for storing the pseudo-absence data 
  pseudoabs.biomod.table <- matrix(FALSE,
                                   nrow = nrow(Fin_occ_var) + runs.PA * number.PA,,
                                   ncol = runs.PA) 
  colnames(pseudoabs.biomod.table) <- paste0("PA", 1:runs.PA) # Column names for each PA run
  pseudoabs.biomod.table[1:nrow(Fin_occ_var), ] <- TRUE # Set the presence points to TRUE
  
  pseudo_only <- data.frame()
  # Loop to sample pseudo-absence points for each run (outside the convex hull and presence pixels)
   for(PA in 1:runs.PA){
    
     # # Sample pseudo-absence points outside the convex hull
     cursp.pseudoabs <- sample(which(!cursp.inhull & !presencepixels), 
                               size = number.PA, # Same number as the presence points
                               replace = FALSE) # No replacement, unique cells
     
     # Pseudo-absence data to the environmental data (NA as values in the "Observed" column)
     cursp.rundata <- rbind.data.frame(
       cursp.rundata,
       data.frame(Observed = NA, # Observed is NA for pseudo-absence points
                  envir.space$unique.conditions.in.env[cursp.pseudoabs, ])) # Environmental conditions for pseudo-absences
     
     pseudo_only <- rbind(pseudo_only,envir.space$coords[cursp.pseudoabs, ])
     # Append the coordinates of pseudo-absence points to the coordinate data
     cursp.xy <- rbind(cursp.xy,envir.space$coords[cursp.pseudoabs, ])
     
     # Update the pseudo-absence table for the current run (mark rows as TRUE for pseudo-absence points)
     pseudoabs.biomod.table[(nrow(Fin_occ_var) + 1 + (PA - 1) * number.PA):
                              (nrow(Fin_occ_var) + PA * number.PA), PA] <- TRUE
   }
  
  
  #### Pseudo absence data for visualisation 
  pseudo_only_vis <- rasterize(as.matrix(pseudo_only), Rastack)
  pseudo_only <- as.data.frame(pseudo_only_vis, xy=T, na.rm = FALSE)

  # Combine occurrences with the baseline raster values
  pseudo_only <- cbind(pseudo_only,Rastab) # Add variable values
  pseudo_only <- pseudo_only[complete.cases(pseudo_only), ] # Remove occurrences outside land
  pseudo_only <- pseudo_only[, c("x","y")]
  
  #### Convexhull Species
  conv_hull <- cbind(data.frame(cursp.inhull), envir.space$coords)
  conv_hull_filtered <- conv_hull[conv_hull$cursp.inhull != FALSE, ]
    
  ### Visualisation data
  plot(Rastack[[1]])
  points(conv_hull_filtered[ , c("x", "y")], pch = 20, cex = 0.5, col = "blue")
  points(pseudo_only[ , c("x", "y")], pch = 20, cex = 0.5, col = "black")
  points(Fin_occ_var[ , c("x", "y")], pch = 20, cex = 0.5, col = "red")
  
  
  # Create the output directory for models if it doesn't exist
  save_dir <- paste0("models/", Sp)
  if(!dir.exists(save_dir)) {
    dir.create(save_dir)
  }
  
  # BIOMOD2 formating data
  run_data <- BIOMOD_FormatingData(
    resp.name = Sp, 
    resp.var = cursp.rundata$Observed, # Observed data for the species (presence/absence values)
    expl.var = cursp.rundata[, -1],   # Predictive variables (species-specific raster stack)
    dir.name = save_dir,  # Directory for saving the models
    resp.xy = cursp.xy,     # XY coordinates of the presence and pseudo-absence points
    PA.strategy = 'user.defined',
    PA.user.table = pseudoabs.biomod.table)
  
  saveRDS(run_data, file = paste0("models/", Sp, "/run_data.RDS"))
}
