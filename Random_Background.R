# ####################### Occurrence filtering & Convexhull
rm(list=ls())
library(terra)
library(openxlsx)
library(raster)
library(biomod2)

# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp


# Change baseline Spat raster as raster
baseline_raster <- as(Rastack, "Raster")

# Transform raster as dataframe
Rastab <- as.data.frame(baseline_raster, xy=T, na.rm = FALSE) # Take one layer of baseline raster

## 1. Occurence filtering
### Remove occurrences coordinates with no correspondence (outside land and false ones)
### Remove occurences with no variable data

i <- 1
# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  Occu <- read.xlsx(paste0("data/raw/occurences/Occurences_", Sp, ".xlsx")) # Load occurrences for this species

  # Proceed only if there are more than 10 occurrences per environmental variable
  if (length(Occu$x) > 60){
    
    
    # Rasterize occurrences to align with baseline raster
    Occu_r <- rasterize(as.matrix(Occu[,1:2]), baseline_raster)
    Occu_r <- as.data.frame(Occu_r, xy=T, na.rm = FALSE) # Convert to data frame with coordinates
    
    Occu_r <- cbind(Occu_r, Rastab) # Add variable values
    Occu_r <- Occu_r[complete.cases(Occu_r), ] # Remove occurrences outside land and with no variable data
    Occu_r <- Occu_r[, c(1,2)]
    Occu_r[3] <- 1 # Assign presence = 1 for each occurrence
    colnames(Occu_r) <- c("x","y","Observed") # Rename columns
    
    
    ##### Random pseudo absences
    
    runs_PA <- 5 # pseudo absence runs
    nb_PA <- 10000 # pseudo absence numbers
    
    # Create the output directory for models if it doesn't exist
    save_dir <- paste0("models/", Sp)
    if(!dir.exists(save_dir)) {
      dir.create(save_dir)
    }
    
    # Pseudo absence generations
    run_data <- BIOMOD_FormatingData(resp.name = Sp, # Species name
                                     resp.var = Occu_r[, 3], # Response variable (Occurrences)
                                     expl.var = Rastack, # Explicatives variables 
                                     dir.name = save_dir, # Folder to save
                                     resp.xy = Occu_r[, c(1,2)], # Occurrences coordinates
                                     PA.nb.rep = runs_PA, # pseudo absence runs
                                     PA.nb.absences = nb_PA,  # pseudo absence numbers
                                     PA.strategy = 'random') # Selection strategy of pseudo absences
    
    saveRDS(run_data, file = paste0("models/", Sp, "/run_data_rdbg.RDS"))
  }
}