library(openxlsx)

# Pseudo abcences generation

Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$x

for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  Occu <- read.xlsx(paste0("data/Filtered_occurences/Occ&Var_", Sp, ".xlsx")) # Load species occurrences

number.PA <- nrow(Occu) # Number of pseudo absence = Number of presence
runs.PA <-5 # Run of pseudo absences

pseudoabs.biomod.table <- matrix(FALSE,
                                 nrow = nrow(Occu) + runs.PA * number.PA,
                                 ncol = runs.PA) # creation table vide qui contient 5 pre run de pseudo-absence : meme nmbre pseudo-absence que de presence
colnames(pseudoabs.biomod.table) <- paste0("PA", 1:runs.PA)
pseudoabs.biomod.table[1:nrow(Occu), ] <- TRUE ## = presence

# presence + conditions VA
cursp.rundata <- Occu[,-c(1,2)] ## Environmental variables
cursp.xy <- Occu[, c("x", "y")] # coordonnées d'occurence des especes

plot = F

}


cursp.xy <- species.occ.filt[, c("x", "y")] # coordonnées d'occurence des especes

plot = F

for(PA in 1:runs.PA){
  
  try({
    # Sampling pseudo-absences outside the convex hull
    # AND not on presence points
    cursp.pseudoabs <- sample(which(!cursp.inhull & !presencepixels), 
                              size = number.PA,
                              replace = FALSE)
    
    cursp.rundata <- rbind.data.frame(cursp.rundata,
                                      data.frame(Observed = NA,
                                                 envir.space$unique.conditions.in.env[cursp.pseudoabs, ]))
    
    cursp.xy <- rbind(cursp.xy,
                      data.frame(xyFromCell(Rastack[[1]],cursp.pseudoabs)))
    
    pseudoabs.biomod.table[(nrow(cur.sp.pixels) + 1 + (PA - 1) * number.PA):
                             (nrow(cur.sp.pixels) + PA * number.PA), PA] <- TRUE
    
  })
  
  if(plot){
    plotconvexhull(allenvpixels = as.data.frame(env.space$unique.conditions.in.env),
                   cursppixels = cur.sp.pixels.filt,
                   pseudoabs = cursp.pseudoabs,
                   curspinhull = cursp.inhull)
  }
}

results <- list(occurrence.environment.matrix = cursp.rundata,
                pseudoabs.biomod.table = pseudoabs.biomod.table, # pour chaque set pseudo abs : occurences + le set de pseudo abs
                xy.coordinates = cursp.xy) ## creation liste avec coordonnées des occurences + pseudo abs (5 car 5 echantillonages) 

filtered.records <- results

coorxy <- filtered.records$xy.coordinates

occurrences <- filtered.records$occurrence.environment.matrix

PATable<-filtered.records$pseudoabs.biomod.table

# Chemin de sauvegarde
save_dir <- paste0(noms_dossiers_bis[j], "/", VecSp[k])
# dir.create(save_dir, recursive = TRUE)  # Création du dossier si nécessaire

# Formatage des données pour BIOMOD2
run_data <- BIOMOD_FormatingData(
  resp.name = VecSp[k], 
  resp.var = occurrences$Observed, 
  expl.var = Rastack,   # Variables prédictives (rasterstack propre à l'espèce)
  dir.name = save_dir,  # Dossier de stockage des modèles
  resp.xy = coorxy,     # Coordonnées xy des présences et pseudo-absences
  PA.strategy = 'user.defined',
  PA.user.table = filtered.records$pseudoabs.biomod.table
)