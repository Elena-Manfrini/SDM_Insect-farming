library(terra)
library(openxlsx)
library(xgboost)
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
# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  
  # Occurrence and environmental values for the species
  Fin_occ_var <- read.xlsx(paste0("data/filtered_occurences/Occ&Var_", Sp, ".xlsx"))
  
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
    test <- cursp.rundata[-calibration, ]
    
    print(paste("Training set size:", dim(train)[1]))
    print(paste("Test set size:", dim(test)[1]))
    
    
    dtrain <- xgb.DMatrix(data = as.matrix(train[,-1]), label = train$Observed)
    dtest <- xgb.DMatrix(data = as.matrix(test[,-1]), label = test$Observed)
    
    # Define parameter grid for tuning
    param_grid <- expand.grid(
      nrounds = c(300,500, 1000),     # Number of boosting rounds
      max_depth = c(5, 7, 9),                  # Maximum depth of a tree
      eta = c(0.1, 0.01, 0.001),               # Learning rate
      gamma = c(0, 1, 5),                      # Minimum loss reduction to make a split
      colsample_bytree = 0.8,       # Subsampling ratio of columns
      min_child_weight = 1,          # Minimum sum of instance weight needed in a child
      subsample = 1
    )
    
    best_auc <- 0
    best_params <- list()
    
    # Hyperparameter tuning using grid search
    for (i in 1:nrow(param_grid)) {
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        max_depth = param_grid$max_depth[i],
        eta = param_grid$eta[i],
        gamma = param_grid$gamma[i],
        colsample_bytree = param_grid$colsample_bytree[i],
        min_child_weight = param_grid$min_child_weight[i],
        subsample = param_grid$subsample[i]
      )
    
      # Cross-validation for AUC
      cv_results <- xgb.cv(
        params = params,
        data = dtrain,
        nrounds = 1000,
        nfold = 5,
        stratified = TRUE,
        verbose = FALSE,
        early_stopping_rounds = 10
      )
      
      mean_auc <- max(cv_results$evaluation_log$test_auc_mean)
      if (mean_auc > best_auc) {
        best_auc <- mean_auc
        best_params <- params
      }
    }
    
    # Train the final model with the best parameters
    final_model <- xgb.train(
      params = best_params,
      data = dtrain,
      nrounds = 1000,
      watchlist = list(train = dtrain, test = dtest),
      verbose = TRUE
    )
    
    # Predict and evaluate
    predictions <- predict(final_model, newdata = dtest)
    pred_labels <- ifelse(predictions > 0.5, 1, 0)
    confusion_matrix <- table(predicted = pred_labels, actual = test$Observed)
    print(confusion_matrix)
    
    auc_score <- auc(roc(test$Observed, predictions))
    print(paste("Best AUC:", auc_score))
}
