rm(list=ls())
library(openxlsx)
library(blockCV)
library(biomod2)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)

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

  # 2.1.a Tune GAM model 
  ## Directly using biomod = too long
  
  # tuned.GAM <- bm_Tuning(model = 'GAM',
  #                        tuning.fun = "gam",
  #                        params.train = list(
  #                          GAM.select = c(TRUE, FALSE),  # Test both with and without variable selection
  #                          GAM.method = c("REML", "GCV.Cp") # ,  # Use REML and GCV.Cp as smoothing selection methods
  #                          # GAM.span = c(0.3, 0.5, 0.7),  # Test different smoothing spans
  #                          # GAM.degree = c(1, 2)  # Use polynomial degrees 1 and 2
  #                        ),
  #                        metric.eval = "ROC",  # Evaluate model performance with ROC
  #                        bm.options = BIOMOD.options.default("GAM", "binary", "mgcv", "gam"),  # Define GAM modeling options
  #                        bm.format = proj_names,  # Pass the formatted projection data
  #                        calib.lines = table_cv  # Use the cross-validation dataset
  #                        )
  
  # After tuning it on a separate Rscript (GAM_param)
  # Define GAM model parameters
  user.GAM <- list('_allData_allRun' = list(select = FALSE, method = "GCV.Cp", 
                                            formula = Hermetia.illucens ~ s(bio5, k=15, sp = 0.01) + s(hurs_min, k=15, sp = 0.01) + s(npp, k=15, sp = 0.01) + s(globalCropland_2010CE,k=15, sp = 0.01)))
  
  # 2.1.b Tune XGBOOST model 
  # After tuning it on a separate Rscript (XGBoost_param)
  # user.XGBOOST <- list('_allData_allRun' = list(nrounds = 500, max.depth = 4))
  
  # 2.1.c Tune RF model
  # After tuning it on a separate Rscript (XGBoost_param)
  # user.RF <- list('_allData_allRun' = list(ntrees = 1000, mtry = 2))
  
  # 2.1.d Combine all tuned models
  user.val <- list(
    # XGBOOST.binary.xgboost.xgboost = user.XGBOOST,
                   # RF.binary.randomForest.randomForest = user.RF,
                   GAM.binary.mgcv.gam = user.GAM)

  # # # Define modeling options
  model_parameters <- bm_ModelingOptions(
    data.type = "binary", # Binary classification (presence/absence)
    models = c(
      # 'RF', 'XGBOOST','MAXNET', 
               'GAM'
               # , 'GLM'
               ),  # Full list of models
    strategy = "user.defined",
    user.val = user.val,
    bm.format = proj_names, # Input model data
    calib.lines = table_cv # Cross-validation table
  )

  ### 2.2 Invividual models with tuned and default model parameters

  # # # Run modeling with specified parameters
  myBiomodModelOut <- BIOMOD_Modeling(
    proj_names,
    modeling.id = "Modeling",
    models = c(
      # 'RF', 'XGBOOST','MAXNET',
               'GAM'
               # , 'GLM'
               ),
    OPT.strategy = "user.defined",
    OPT.user = model_parameters,
    CV.strategy = 'user.defined',
    CV.user.table = table_cv,
    CV.do.full.models = FALSE,
    var.import = 10, # Number of variable importance metrics to calculate
    metric.eval = c('BOYCE'), # Evaluation metric (Boyce index)
    do.progress = TRUE)

  
  # # # Run modeling with specified parameters
  # myBiomodModelOut_default <- BIOMOD_Modeling(
  #   proj_names,
  #   modeling.id = "Modeling",
  #   models = c(
  # # 'RF', 'XGBOOST','MAXNET',
  # 'GAM'),
  #   OPT.strategy = "default", # Try to optimize compare to default
  #   CV.strategy = 'user.defined',
  #   CV.user.table = table_cv,
  #   CV.do.full.models = FALSE,
  #   var.import = 10, # Number of variable importance metrics to calculate
  #   metric.eval = c('BOYCE'), # Evaluation metric (Boyce index)
  #   do.progress = TRUE)
  
  # Save the modeling output
  saveRDS(myBiomodModelOut, file = paste0("models/", Sp, "/output_models_GAM.rds"))
  
  
  ## 2.2.a Response curves
  
  # Load pre-trained model outputs for the species
  myBiomodModelOut <- readRDS(paste0("models/", Sp, "/output_models_GAM.RDS"))
  
  # Get the evaluation results for each model
  evals <- get_evaluations(myBiomodModelOut)
  
  # Filter models with a Boyce index greater than 0.7
  Choosen_Model <- evals %>% filter(validation > 0.7)
  
  # Extract the names of selected models
  selected_models <- Choosen_Model$full.name # Get the full names of selected models
  
  # Save the selected models
  saveRDS(selected_models, file = paste0("models/", Sp, "/selected_models_GAM.rds"))
  
  # Plot response curves for the selected models
  resp <- bm_PlotResponseCurves(bm.out = myBiomodModelOut,
                              models.chosen = selected_models,
                              # show.variables = cur_vars,http://127.0.0.1:37987/graphics/plot_zoom_png?width=1029&height=809
                              fixed.var = "mean", # Fix the variable at the mean value
                              data_species = proj_names@data.species,
                              do.plot = FALSE, # Do not plot here (only generate the response table)
                              do.progress = FALSE)$tab
  
  colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")
  
  # Extract the model type from the 'Model' column
  resp <- resp %>%
    mutate(Model_Type = sub(".*_(GLM|GAM|MAXNET|RF|XGBOOST|MAXNET)$", "\\1", Model))
  
  
  response <- ggplot(resp, aes(x = Var.value, y = Response)) + 
    geom_line(alpha = 0.2, aes(group = Model, color = Model_Type)) +  # Use Model_Type for color
    stat_smooth(color = "black") + #Trend line
    facet_wrap(~Variable, scales = "free_x") + 
    theme_bw() +  # Black-and-white theme
    ylim(0, 1.1) + 
    xlab("Variable value") +
    scale_color_manual(values = c(
      # 'RF' = '#26c031', 
      # 'XGBOOST' = '#38599f', 
      # 'MAXNET' = '#db1010',
      # 'GLM' = 'lightblue',
      'GAM' = 'red'
      )) +  # Custom color scale
    labs(color = "Model")

   
  # Create the output directory for figures if it doesn't exist
  save_dir <- paste0("output/Relation_Environnement")
  if(!dir.exists(save_dir)) {
  dir.create(save_dir)
  }
  
  # Save the plot 
  ggsave(str_c(save_dir, "/",Sp,"_GAM.jpeg",sep= ""),
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
    scale_color_manual(values = c(
      'RF' = '#26c031', 
      'XGBOOST' = '#38599f', 
      'MAXNET' = '#db1010')) +  # Custom color scale
    geom_jitter (alpha = 0.2, aes(col = Model))
  
  # Create the output directory for figures if it doesn't exist
  save_dir <- paste0("output/Variable_Importance")
  if(!dir.exists(save_dir)) {
    dir.create(save_dir)
  }
  
  # Save the plot 
  ggsave(str_c(save_dir, "/",Sp,"_hfp.jpeg",sep=""),
         var_imp,
         dpi = 2000,
         bg = NULL,
         width = 11,
         height = 8.5,
         units = "in")
  
}
