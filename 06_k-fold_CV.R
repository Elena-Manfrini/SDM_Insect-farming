library(openxlsx)
library(blockCV)
library(biomod2)
library(dplyr)
library(tidyverse)

library(ecospat)

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

i <- 1
# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Get the current species name
  
  # Read model projection data for the species
  proj_names <- readRDS(paste0("models/", Sp, "/run_data.RDS"))
  
  ############# 1. K-fold crossvalidation
  
  # Perform K-fold cross-validation with 5 folds
  table_cv <- bm_CrossValidation(
    bm.format = proj_names, # Model projections to validate
    strategy = "kfold", # Cross-validation strategy
    k = 5, # Number of folds for cross-validation
    nb.rep = 1, # Number of repetitions
    balance = "both") # Option to balance presence and absence data
  
  # Extract the calibration summary from the cross-validation results
  calib_summary <-
    summary(proj_names, calib.lines =  table_cv) %>%
    filter(dataset == "calibration") # Filter calibration dataset for the cross-validation
  
  ############# 2. Model formating
  
  ### 2.1 Model preparation
  
  # Create empty lists to store parameter configurations for different models
  RF_param_list <- list()
  XGBOOST_param_list <- list()
  MAXNET_param_list <- list()
  
  # Loop over each cross-validation run to adjust model parameters
  for (cvrun in 1:nrow(calib_summary)) {
    prNum <- calib_summary$Presences[cvrun] # Number of presence points
    bgNum <- calib_summary$Pseudo_Absences[cvrun] # Number of pseudo abscence points
    
    # Random Forest model parameters
    RF_param_list[[paste0("_", calib_summary$PA[[cvrun]], "_", calib_summary$run[[cvrun]])]] <- list(
       ntree = 2000,  # Number of trees in the forest
       mtry = floor(sqrt(ncol(proj_names@data.env.var))), # number of variables per tree
       sampsize = c("0" = bgNum, "1" = bgNum), # Balanced sampling of presence and background
       replace = TRUE
     )
     
    # XGBOOST model parameters
    # Adjust weights for the XGBOOST model based on presence and background points
    wt <- ifelse(proj_names@data.species == 1, 1, prNum / bgNum) #### toujours 1 : meme poids pour les pseudo abs que pour les presences
    
  XGBOOST_param_list[[paste0("_", calib_summary$PA[[cvrun]], "_", calib_summary$run[[cvrun]])]] <- list(
    nrounds = 1000,  # number of iterations
    eta = 0.1,  # learning rate
    max_depth = 7,  # depth of trees
    subsample = 1,  # Use 90% of the data for each tree
    objective = "binary:logistic", # Binary logistic regression
    gamma = 0,  # Regularization to avoid overfitting
    colsample_bytree = 0.8,  # Fraction of features to sample per tree
    min_child_weight = 1, # Minimum sum of instance weight
    weight = wt, # Use calculated weights for each data point
    verbose = 0
  )
  
  # # MAXNET model parameters
   MAXNET_param_list[[paste0("_", calib_summary$PA[[cvrun]], "_", calib_summary$run[[cvrun]])]] <- list(
     l1_regularizer = 0.1,  # L1 regularization to prevent overfitting
     l2_regularizer = 0.1,  # L2 regularization to prevent overfitting
     use_sgd = TRUE, # Use stochastic gradient descent
     set_heldout = 0, # No held-out data for validation
     verbose = TRUE
   )
  }
  
  # Define modeling options
  model_parameters <- bm_ModelingOptions(
    data.type = "binary", # Binary classification (presence/absence)
    models = c('RF', 'XGBOOST','MAXNET', 'GAM', 'GLM', 'GBM'),  # Full list of models
    strategy = "user.defined",
    user.base = "default",
    user.val = list(
      RF.binary.randomForest.randomForest = RF_param_list,
      XGBOOST.binary.xgboost.xgboost = XGBOOST_param_list,
      MAXNET.binary.maxnet.maxnet = MAXNET_param_list), 
    bm.format = proj_names, # Input model data
    calib.lines = table_cv # Cross-validation table
  ) 
  
  ### 2.2 Invividual models
  
  # Run modeling with specified parameters
  myBiomodModelOut <- BIOMOD_Modeling(
    proj_names,
    modeling.id = "Modeling",
    models = c('RF', 'XGBOOST','MAXNET', 'GAM', 'GLM', 'GBM'),
    OPT.strategy = "user.defined",
    OPT.user = model_parameters,
    CV.strategy = 'user.defined',
    CV.user.table = table_cv,
    CV.do.full.models = FALSE,
    var.import = 10, # Number of variable importance metrics to calculate
    metric.eval = c('BOYCE'), # Evaluation metric (Boyce index)
    do.progress = TRUE)
  
  # Save the modeling output
  saveRDS(myBiomodModelOut, file = paste0("models/", Sp, "/model_output.rds"))
  
  
  ## 2.2.a Reponse curves
  
  # Get the evaluation results for each model
  evals <- get_evaluations(myBiomodModelOut)
  
  # Filter models with a Boyce index greater than 0.7
  Choosen_Model <- evals %>% filter(validation > 0.7)
  
  # Extract the names of selected models
  selected_models <- Choosen_Model$full.name # Get the full names of selected models
  
  # Save the selected models
  saveRDS(selected_models, file = paste0("models/", Sp, "/selected_models.rds"))
  
  # Plot response curves for the selected models
  resp <- bm_PlotResponseCurves(bm.out = myBiomodModelOut,
                              models.chosen = selected_models,
                              # show.variables = cur_vars,http://127.0.0.1:37987/graphics/plot_zoom_png?width=1029&height=809
                              fixed.var = "mean", # Fix the variable at the mean value
                              data_species = proj_names@data.species,
                              do.plot = FALSE, # Do not plot here (only generate the response table)
                              do.progress = FALSE)$tab
  
  colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")
  
  # Plot the response curves for each model
  response <- ggplot() +
    geom_line(data = resp[grep('RF', resp$Model), ],
              aes(x = Var.value, y = Response, color = 'RF'), alpha = 0.2) +
    geom_line(data = resp[grep('XGBOOST', resp$Model), ],
              aes(x = Var.value, y = Response, color = 'XGBOOST'), alpha = 0.2) +
    geom_line(data = resp[grep('MAXNET', resp$Model), ], 
              aes(x = Var.value, y = Response, color = 'MAXNET'), alpha = 0.2) +
    geom_line(data = resp[grep('GAM', resp$Model), ], 
              aes(x = Var.value, y = Response, color = 'GAM'), alpha = 0.2) +
    geom_line(data = resp[grep('GLM', resp$Model), ], 
              aes(x = Var.value, y = Response, color = 'GLM'), alpha = 0.2) +
    geom_line(data = resp[grep('GBM', resp$Model), ], 
              aes(x = Var.value, y = Response, color = 'GBM'), alpha = 0.2) +
    facet_wrap(~ Variable, scales = "free_x") +
    theme_bw() +
    ylim(0, 1.1) +
    xlab("Variable value") +
    scale_color_manual(values = c('RF' = 'darkgreen', 
                                  'XGBOOST' = 'darkblue', 
                                  'MAXNET' = 'darkred', 
                                  'GAM' = 'purple', 
                                  'GLM' = 'orange', 
                                  'GBM' = 'cyan')) +  # Custom color scale
    labs(color = "Model")
  
  # Create the output directory for figures if it doesn't exist
  save_dir <- paste0("output/Relation_Environnement")
  if(!dir.exists(save_dir)) {
  dir.create(save_dir)
  }
  
  # Save the plot 
  ggsave(str_c(save_dir, "/",Sp,".jpeg",sep=""),
         response,
       dpi = 2000,
       bg = NULL,
       width = 11,
       height = 8.5,
       units = "in")
  
  
  ## 2.2.b Variable importance
  
  gg_varimp <- get_variables_importance(myBiomodModelOut) # techniques de model pour le run + Va + importance calculée des VA
  
  
  colnames(gg_varimp) <- c("id", 
                           "PA.Run",
                           "CV.Run", "Model", 
                           "Variable", 
                           "VI.run",
                           "Variable.importance")
  
  var_imp <- ggplot(gg_varimp, aes(y = Variable, x = Variable.importance)) + # ordonne le graph en fonction de celle 
    #qui a la plus forte mediane : utilise les niveaux de facteur 
    geom_boxplot() + 
    theme_bw() +
    ggtitle(Sp) + 
    scale_color_manual(values = c('RF' = 'darkgreen', 
                                  'XGBOOST' = 'darkblue', 
                                  'MAXNET' = 'darkred', 
                                  'GAM' = 'purple', 
                                  'GLM' = 'orange', 
                                  'GBM' = 'cyan')) +  # Custom color scale
    geom_jitter (alpha = 0.2, aes(col = Model))
  
  # Create the output directory for figures if it doesn't exist
  save_dir <- paste0("output/Variable_Importance")
  if(!dir.exists(save_dir)) {
    dir.create(save_dir)
  }
  
  # Save the plot 
  ggsave(str_c(save_dir, "/",Sp,".jpeg",sep=""),
         var_imp,
         dpi = 2000,
         bg = NULL,
         width = 11,
         height = 8.5,
         units = "in")
  
}
