library(terra)
library(openxlsx)
library(biomod2)

# Load Environmental Space
envir.space <- readRDS("data/Environmental_Space.rds")

# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  
  # Coordinates & Environemental values of species
  Fin_occ_var <- read.xlsx(paste0("data/filtered_occurences/Occ&Var_", Sp, ".xlsx"))
  
  # ConvexHull
  cursp.inhull <- readRDS(paste0("data/convexhull/", Sp, "_cursp.inhull.rds"))
  # Occurences inside ConvexHull
  presencepixels <- readRDS(paste0("data/convexhull/", Sp, "_presencepixels.rds"))
  
  number.PA <- nrow(Fin_occ_var)
  runs.PA <-5
  
  cursp.rundata <- Fin_occ_var[,-c(1:2)] # Environmental conditions
  
  pseudoabs.biomod.table <- matrix(FALSE,
                                   nrow = nrow(Fin_occ_var) + runs.PA * number.PA,,
                                   ncol = runs.PA) # creation table vide qui contient 5 pre run de pseudo-absence : meme nmbre pseudo-absence que de presence
  colnames(pseudoabs.biomod.table) <- paste0("PA", 1:runs.PA)
  pseudoabs.biomod.table[1:nrow(Fin_occ_var), ] <- TRUE
  
  # Occurrences coordinates
  cursp.xy <- Fin_occ_var[, c("x", "y")] 
  
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
    
    pseudoabs.biomod.table[(nrow(Fin_occ_var) + 1 + (PA - 1) * number.PA):
                             (nrow(Fin_occ_var) + PA * number.PA), PA] <- TRUE
  }

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
    resp.var = cursp.rundata$Observed, 
    expl.var = cursp.rundata[, -1],   # Variables prédictives (rasterstack propre à l'espèce)
    dir.name = save_dir,  # Dossier de stockage des modèles
    resp.xy = cursp.xy,     # Coordonnées xy des présences et pseudo-absences
    PA.strategy = 'user.defined',
    PA.user.table = pseudoabs.biomod.table)
  
  saveRDS(run_data, file = paste0("models/", Sp, "/run_data.RDS"))
}