######################## Occurrence filtering & Convexhull ########################
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

out_in_ConvHull_All <- data.frame(Species = character(),
                              NbOcc_InConv = numeric(),
                              NbOcc_OutConv = numeric(),
                              stringsAsFactors = FALSE)


# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  Occu <- read.xlsx(paste0("data/raw/occurences/Occurences_", Sp, ".xlsx")) # Load occurrences for this species
 
   ############# 1. Occurrences filtering
  
    ### 1.1 Remove occurrences outside land
    
    # Proceed only if there are more than 10 occurrences per environmental variable
    if (length(Occu$x) > 50){
    # Rasterize occurrences to align with baseline raster
    Occu_r<- rasterize(as.matrix(Occu[,1:2]), baseline_raster)
    Occu_r <- as.data.frame(Occu_r, xy=T, na.rm = FALSE) # Convert to data frame with coordinates
    
    # Combine occurrences with the baseline raster values
    Occu_r <- cbind(Occu_r,Rastab) # Add variable values
    Occu_r <- Occu_r[complete.cases(Occu_r), ] # Remove occurrences outside land
    Occu_r <- Occu_r[, c("x", "y", "layer")]
    Occu_r[3] <- 1 # Assign presence = 1 for each occurrence
    colnames(Occu_r) <- c("x","y","Observed") # Rename columns
    
    ### 1.2 Remove duplicated occurrences in Environmental Space
    
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
    

    # Final dataframe of occurences and their respective variable values
    Fin_occ_var <- cbind(Occu_2, var.occ_2)
    Fin_occ_var <- Fin_occ_var %>%
      dplyr::select(-ID)

    # Save occurrences to an Excel file
    if(!dir.exists("data/filtered_occurences")) {
      dir.create("data/filtered_occurences,")
    }
    xlsx::write.xlsx(Fin_occ_var, paste0("data/filtered_occurences/Occ&Var_final", Sp, ".xlsx"), row.names = F)


    ############# 2. Convex Hull

    ### 2.1 Convex Hull preparation
    # Replace interval values by their mid points.
    cur.sp.pixels <- sapply(colnames(env.itvl_2), function(x, int.to.replace, replacing.values){
      replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
    }, int.to.replace = env.itvl_2, replacing.values = envir.space$detailed.intervals) ### valeurs variables correspondant aux occurences uniques

    # Remove occurrences outliers for the convex hull calculations (1% per variable)
    convex.hull.interval = 0.01
    

    # Select for each variable 0.5% of the lower extreme and 0.5% of the higher extreme
    outs <- lapply(colnames(var.occ_2), function(x, df)
    {
      qt <- quantile(df[, x], probs = c(0 + convex.hull.interval / 2,
                                        1 - convex.hull.interval / 2))
      return(which(df[, x] <= qt[1] | df[, x] >= qt[2]))
    }, df = var.occ_2)

    outs <- unique(unlist(outs))

    # Remove occurrences outliers from midpoint dataframe
    cur.sp.pixels.filt <- cur.sp.pixels[-outs,]


    ### 2.2 Convex Hull creation

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
    
    ### 2.3 Presence pixel 
    
    # Check where species occurences are located into the environmental space (using midpoints variable values)
    presencepixels <- apply(envir.space$unique.conditions.in.env, 1, paste, collapse = " ")   %in%
      ## For each row, the env values are concatenated into a single string with spaces between them.
      # check if elements of one vector are present in another.
      apply(cur.sp.pixels, 1, paste, collapse = " ")
    #It returns a logical vector (TRUE or FALSE), indicating whether each element in the first vector is found in the second vector

    # Save presence pixels
    saveRDS(presencepixels, paste0("data/convexhull/", Sp, "_presencepixels.rds"))

    } else {
      # Stop the script if the number of rows in Occu is not greater than 60
      stop("The number of occurrences is less than or equal to 50. Exiting the script.")
    }
}
