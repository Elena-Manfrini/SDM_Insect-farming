setwd("~/Documents/Documents/Eléna/Thèse_Eléna/2024-2025/04_SDMs/Reproductibility")

library(terra)
library(sf)
library(Rarity) # plot des correlations multiples
library(virtualspecies)
library(sp)
library(raster)
library(geodata)
library(openxlsx)

###################### pour charger fonction computeEnvCombinations. ######################
computeEnvCombinations <- function(env.stack,
                                   var.intervals,
                                   plot = TRUE,
                                   vars.to.plot = 1:3)
{
  combinations <- as.data.frame(getValues(Rastack[[c(names(Rastack))]]))
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
  
  possible.combs <- lapply(var.intervals, function(x)
    data.frame(interval = cut(x, x, right = FALSE, include.lowest = TRUE),
               mids = midpoints(cut(x, x, right = FALSE, include.lowest = TRUE))))
  
  all.env.pixels <- sapply(colnames(comb.cat), function(x, int.to.replace, replacing.values)
  {
    replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
  }, int.to.replace = comb.cat, replacing.values = possible.combs)
  
  duplicated.cells <- which(duplicated(all.env.pixels))
  
  if(length(duplicated.cells)){
    
    all.env.pixels <- all.env.pixels[-duplicated.cells, ]
    
    all.xy <- all.xy[-duplicated.cells, ]
  }
  
  if(plot)
  {
    plot3d(all.env.pixels[, vars.to.plot])
  }
  return(list(detailed.intervals = possible.combs,
              unique.conditions.in.env = all.env.pixels,
              coords = all.xy))
}


#### package midpoints 

midpoints <- function(x, dp=2){
  lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
  upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
  return(round(lower+(upper-lower)/2, dp))
}
#########################################################################################


############# 1. Variables environnementales  

# Chargement des données environnementales
baseline <- rast("data/donnees_brutes/bioclim/final_baseline.tif")
# Renommer manuellement
names(baseline) <- c("CHELSA_bio5", "CHELSA_bio7", "CHELSA_hurs_min", "CHELSA_hurs_range", "CHELSA_npp", "globalCropland_2010CE")

# Calculer les valeurs environnementales (quantiles)
values <- values(baseline) # recuperation valeurs ; reprend les NA des oceans
combinations <- as.data.frame(values)
combinations <- combinations[complete.cases(combinations), ] #numero de pixel

# Calcul des quantiles pour chaque variable : separation des valeurs en 100 : grille pour localiser les occurences
Bio5 <- quantile(combinations[,1], probs=seq(0,1,0.01))
Bio7 <- quantile(combinations[,2], probs=seq(0,1,0.01))
Hurs_min <- quantile(combinations[,3], probs=seq(0,1,0.01))
Hurs_range <- quantile(combinations[,4], probs=seq(0,1,0.01))
Npp <- quantile(combinations[,5], probs=seq(0,1,0.01))
croplands <- quantile(combinations[,6], probs=seq(0,1,0.01))

# Créer une liste d'intervalles environnementaux uniques pour chaque variable : Construction du quadrillage en 5 dimensions
intervals <- list(
  Bio5 = sort(unique(as.numeric(Bio5))), # unique = une occurence par carré de la grille
  Bio7 = sort(unique(as.numeric(Bio7))),
  Hurs_min = sort(unique(as.numeric(Hurs_min))),
  Hurs_range = sort(unique(as.numeric(Hurs_range))),
  Npp = sort(unique(as.numeric(Npp))),
  croplands = sort(unique(as.numeric(croplands))) #### Uniquement 32 au lieu de 101 valeurs
)

# Assigner des noms aux intervalles en fonction des variables du stack
names(intervals) <- names(baseline)

#intervals$globalCropland_2010CE <- seq(min(intervals$globalCropland_2010CE),
                                       #max(intervals$globalCropland_2010CE), 
                                       #length.out = 101)

# Calcul des combinaisons environnementales : creation espace environnemntale
envir.space <- computeEnvCombinations(
  env.stack = baseline[[c(names(baseline))]],
  var.intervals = intervals, #grille definie
  plot = F,
  vars.to.plot = 1:6) ##### Erreur du à Cropland qui ne possède que 32 intervals ?

############# 2. Ajout occurrences

# Recupération des occurences des espèces 
# Recuperer toutes les variables
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$x
baseline_raster <- as(baseline, "Raster")


for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]]
  Occu <- read.xlsx(paste0("./data/donnees_brutes/occurences/Occurences_", Sp, ".xlsx"))
  
  # Si on obtient plus de 50 occurrences, on continue
  if (length(Occu$x) > 50){
    
    # Rasteriser les occurrences
    Occu <- rasterize(Occu[,1:2],baseline_raster)
    Occu <- as.data.frame(Occu, xy=T)
    
    # Combiner avec le stack climatique : avoir un raster final avec 0 NA
    Rastab <- as.data.frame(baseline_raster[[1]],xy=T)
    Occu <- cbind(Occu,Rastab[[3]])
    Occu <- Occu[complete.cases(Occu[ , c(colnames(Occu))]), ] ## supprimer toutes les occurences qui sont situées hors du sol (mer)
    Occu <- Occu[,-c(4)]
    Occu$layer <- rep(1,length(Occu$x))
    colnames(Occu) <- c("x","y","Observed")
  }
  
  var.intervals <- intervals
  
  species.occurrence <- Occu
  
  cur.occ <- raster::extract(baseline_raster[[c(names(baseline_raster))]],
                             species.occurrence[, c("x", "y")]) #### extraire les valeurs des VA env aux points d'occurences
  
  env.comb.cursp <- sapply(colnames(cur.occ), function(x, combs, seqs.)
    cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
    combs = cur.occ, seqs. = var.intervals) #### Dans quelles cases vont se fixer chacune des occurences par rapport aux intervals fournis initiallement ? Où vont se loger les occurences dans l'espace env
  
  duplicated.cells <- which(duplicated(env.comb.cursp)) # cellules où occurences sont dupliquées
  
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
  
  if(length(duplicated.cells)>0){
    species.occ.filt <- species.occurrence[-duplicated.cells, ] ## on enleve les cellules dupliquées
    env.comb.cursp <- env.comb.cursp[-duplicated.cells, ]
  }else{
    species.occ.filt<-species.occurrence 
  } 
  
  cur.sp.pixels <- sapply(colnames(env.comb.cursp), function(x, int.to.replace, replacing.values){
    replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
  }, int.to.replace = env.comb.cursp, replacing.values = envir.space$detailed.intervals) ### valeurs variables correspondant aux occurences uniques
  
  if(any(species.occurrence$Observed == 0)){
    
    warning("Your occurrence dataset has absences and you are bout to generate pseudoabsences; make sure it
              is correct (or disable generate.PA)")
  }
  
  number.PA = "Number of presences"
  
  if(number.PA != "Number of presences"){
    
    if(!is.numeric(number.PA)){
      stop("number.PA must either be a numeric value or 'Number of presences'")
    }
    
  }else{
    number.PA <- length(which(species.occ.filt$Observed == 1)) ## Nbre abscence = nbre presence
    number.PA
  }
  
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
  
  
}



