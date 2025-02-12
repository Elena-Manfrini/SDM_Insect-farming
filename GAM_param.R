rm(list=ls())
library(terra)
library(openxlsx)
library(xgboost)
library(caret)
library(pROC)
library(ggplot2)
library(viridis)

# Load Environmental Space
envir.space <- readRDS("data/Environmental_Space.rds")

# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

i <-  1
# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  
  # Occurrence and environmental values for the species
  Fin_occ_var <- read.xlsx(paste0("data/filtered_occurences/Occ&Var_final_15_bio6", Sp, ".xlsx"))
  
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
  
  ### We can only tune select and method parameters with biomod2 package
  k_values <- c(10, 15, 25, 50)
  
  # select_values = c(TRUE, FALSE)
  methods <- c("GCV.Cp", "GACV.Cp", "REML", "P-REML", "ML", "P-ML")
  # Create a dataframe to store results
  results <- data.frame(K = integer(), METHOD = character(), AIC = numeric(), Train_AUC = numeric(), Test_AUC = numeric())
  
  for (method in methods) {
    for (k in k_values) {
      # Fit GAM model
      model <- gam(Observed ~ s(bio5, k = k) + s(hurs_min, k = k) + s(bio6, k = best_k) +
                     s(npp, k = k) + s(globalCropland_2010CE, k = k),
                   family = binomial, data = train, method = method, select = FALSE)
      
      # Compute AIC
      model_aic <- AIC(model)
      
      # Predict on training set
      train <- train %>% mutate(predicted_train = predict(model, newdata = train, type = "response"))
      train_roc <- roc(as.numeric(train$Observed), train$predicted_train)
      train_auc <- auc(train_roc)
      
      # Predict on test set
      test <- test %>% mutate(predicted_test = predict(model, newdata = test, type = "response"))
      test_roc <- roc(as.numeric(test$Observed), test$predicted_test)
      test_auc <- auc(test_roc)
      
      # Store results
      results <- rbind(results, data.frame(K = k, METHOD = method, AIC = model_aic, Train_AUC = train_auc, Test_AUC = test_auc))
    }
  }
  
  # Convert k to factor for better visualization
  results$K <- as.factor(results$K)
  
  ggplot(results, aes(x = K)) +
    geom_line(aes(y = Train_AUC, group = METHOD, color = "Train AUC"), size = 1) +
    geom_line(aes(y = Test_AUC, group = METHOD, color = "Test AUC"), size = 1, linetype = "dashed") +
    facet_wrap(~ METHOD) +
    scale_color_manual(values = c("Train AUC" = "blue", "Test AUC" = "red")) +
    labs(title = "Overfitting Check: Training vs. Test AUC",
         x = "Number of Basis Functions (k)",
         y = "AUC",
         color = "Dataset") +
    theme_minimal()
  
  ##### Best model parametrisation 
  
  # Find the best model based on AIC and AUC
  best_model <- results %>% arrange(AIC, desc(Train_AUC)) %>% slice(1)
  best_k <- as.factor(best_model$K)
  best_k <- as.numeric(as.character(best_k))
  
  best_method <- best_model$METHOD
  
  train_f <- train[,-c(7)]
  
  # Train the best model
  best_model_run <- gam(Observed ~ s(bio5, k = best_k) + s(hurs_min, k = best_k) + s(bio6, k = best_k)+
                          s(npp, k = best_k) + s(globalCropland_2010CE, k = best_k),
                        family = binomial, data = train_f, method = best_method, select = FALSE)
  
  par(mfrow = c(2, 2))
  gam.check(best_model_run)
  
  k.check(best_model_run)
  
  
  # Predict on test set using the best model
  test <- test %>% mutate(best_predicted = predict(best_model_run, newdata = test, type = "response"))
  
  # Compute AUC for the best model
  best_roc_curve <- roc(as.numeric(test$Observed), test$best_predicted)
  best_auc <- auc(best_roc_curve)
  