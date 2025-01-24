rm(list=ls())
library(terra)
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(biomod2)
library(dplyr)


# var.intervals <- intervals

###################### Define computeEnvCombinations function ######################
computeEnvCombinations <- function(env.stack,
                                   var.intervals,
                                   plot = TRUE,
                                   vars.to.plot = 1:4)
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

# intervals <- list(
#   bio5 <- seq(floor(min(combinations[, 1])), ceiling(max(combinations[, 1])), length.out = 5),
#   hurs_min <- seq(floor(min(combinations[, 2])), ceiling(max(combinations[, 2])), by = 2),
#   # hurs_mean <- seq(floor(min(combinations[, 3])), ceiling(max(combinations[, 3])), by = 1),
#   # hurs_range <- seq(floor(min(combinations[, 4])), ceiling(max(combinations[, 4])), by = 1),
#   npp <- seq(floor(min(combinations[, 3])), ceiling(max(combinations[, 3])),by = 100)
#   # croplands <- seq(floor(min(combinations[, 6])), ceiling(max(combinations[, 6])), by = 0.01)
#   )

intervals <- list(
  bio5 = seq(min(combinations[, 1]), max(combinations[, 1]), length.out = 80),
  hurs_min = seq(min(combinations[, 2]), max(combinations[, 2]), length.out = 80),
  npp = seq(min(combinations[, 3]), max(combinations[, 3]), length.out = 80),
  # Human_pop_2000 = seq(min(combinations[, 4]), max(combinations[, 4]), length.out = 80)
   # ,
  globalCropland_2010CE = seq(min(combinations[, 4]), max(combinations[, 4]), length.out = 80)
)
  
names(intervals) <- names(Rastack)
  
# Optionally save intervals for each step size as a separate file
saveRDS(intervals, file = paste0("data/intervals.rds"))
  
envir.space <- computeEnvCombinations(
      env.stack = Rastack,
      var.intervals = intervals,
      plot = TRUE,
      vars.to.plot = 1:4
    )
    
# Save environmental space for each step size
saveRDS(envir.space, file = paste0("data/Environmental_Space.rds"))
