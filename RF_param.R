library(terra)
library(openxlsx)
library(randomForest)
library(caret)
library(pROC)

# Load Environmental Space
envir.space <- readRDS("data/Environmental_Space.rds")

# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

best_model <- data.frame()

i <-  1
# Create an empty list to store the rasters for each variable
species.estimation <- list()

# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  
  # Occurrence and environmental values for the species
  Fin_occ_var <- read.xlsx(paste0("data/filtered_occurences/final_Occ&Var_", Sp, ".xlsx"))
  
  # Convex hull and presence pixels data
  cursp.inhull <- readRDS(paste0("data/convexhull/", Sp, "_cursp.inhull.rds"))
  presencepixels <- readRDS(paste0("data/convexhull/", Sp, "_presencepixels.rds"))
  
  
  # Set the number of presence-absence (PA) points and the number of random runs for PA sampling
  number.PA <- nrow(Fin_occ_var) # Number of presence points = Number of pseudo absence points
  runs.PA <- 5 # Number of random runs for pseudo-absence sampling
  
  # Prepare environmental data for the species 
  cursp.rundata <- Fin_occ_var[,-c(1:2)] # Environmental conditions
  
  # Table for storing the pseudo-absence data 
  pseudoabs.biomod.table <- matrix(FALSE,
                                   nrow = nrow(Fin_occ_var) + runs.PA * number.PA,,
                                   ncol = runs.PA) 
  colnames(pseudoabs.biomod.table) <- paste0("PA", 1:runs.PA) # Column names for each PA run
  pseudoabs.biomod.table[1:nrow(Fin_occ_var), ] <- TRUE # Set the presence points to TRUE
  
  # XY coordinates of the species' presence points
  cursp.xy <- Fin_occ_var[, c("x", "y")] 
  
  # Loop to sample pseudo-absence points for each run (outside the convex hull and presence pixels)
  for (PA in 1:runs.PA){
    
    # # Sample pseudo-absence points outside the convex hull
    cursp.pseudoabs <- sample(which(!cursp.inhull & !presencepixels), 
                              size = number.PA, # Same number as the presence points
                              replace = FALSE) # No replacement, unique cells
    
    # Pseudo-absence data to the environmental data (NA as values in the "Observed" column)
    cursp.rundata <- rbind.data.frame(
      cursp.rundata,
      data.frame(Observed = NA, # Observed is NA for pseudo-absence points
                 envir.space$unique.conditions.in.env[cursp.pseudoabs, ])) # Environmental conditions for pseudo-absences
    
    # Append the coordinates of pseudo-absence points to the coordinate data
    cursp.xy <- rbind(cursp.xy,
                      data.frame(xyFromCell(Rastack[[1]], 
                                            cursp.pseudoabs)))
    
    # Update the pseudo-absence table for the current run (mark rows as TRUE for pseudo-absence points)
    pseudoabs.biomod.table[(nrow(Fin_occ_var) + 1 + (PA - 1) * number.PA):
                             (nrow(Fin_occ_var) + PA * number.PA), PA] <- TRUE
    
  }
  
  cursp.rundata$Observed[is.na(cursp.rundata$Observed)] <- 0
  cursp.rundata$Observed <- as.character(cursp.rundata$Observed)
  
  calibration <- sample(nrow(cursp.rundata), 0.7 * nrow(cursp.rundata))
  
  train <- cursp.rundata[calibration, ]
  train$Observed <- as.factor(train$Observed)
  
  test <- cursp.rundata[-calibration, ]
  test$Observed <- as.factor(test$Observed)
  
  print(paste("Training set size:", dim(train)[1]))
  print(paste("Test set size:", dim(test)[1]))
 
    trees <- c(500, 1000, 2000, 3000)
    mtry <- c(1,2)
    
    # Create an empty list to store performance for each combination of trees and mtry
    true.estimation <- list()
    
    # Loop over trees and mtry values
    for (j in 1:length(trees)) {
      ntree <- trees[j]
      
      for (k in 1:length(mtry)) {
        nmtry <- mtry[k]
        
        # Train Random Forest model
        rf_model <- randomForest(
          Observed ~ .,  # Predict 'Observed' based on other variables
          data = train,   # Training dataset
          importance = TRUE,   # Compute variable importance
          mtry = nmtry,    # Number of variables to try at each split
          ntree = ntree    # Number of trees
        )
        
        # Make predictions on the test data
        pred <- predict(rf_model, test)
        
        # Combine actual and predicted values
        class_rf <- as.data.frame(cbind(test$Observed, pred))
        colnames(class_rf)[1] <- "actual_values"
        colnames(class_rf)[2] <- "predicted_values"
        
        # Calculate the accuracy for this combination of ntree and mtry
        accuracy <- (nrow(subset(class_rf, actual_values == 2 & predicted_values == 2)) + 
                       nrow(subset(class_rf, actual_values == 1 & predicted_values == 1))) / 
          nrow(class_rf) * 100
        
        # Store the result with a unique key combining ntree and mtry
        true.estimation[[paste("ntree_", ntree, "_mtry_", nmtry, sep = "")]] <- accuracy
      }
    }
    species.estimation[[Sp]] <- true.estimation
}