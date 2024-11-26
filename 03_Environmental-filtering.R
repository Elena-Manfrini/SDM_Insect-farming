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
}