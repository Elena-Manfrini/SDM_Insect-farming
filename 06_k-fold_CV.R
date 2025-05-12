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

RF_param <- read.xlsx("data/best_RF_param.xlsx")
XGBOOST_param <- read.xlsx("data/best_XGBOOSt_param.xlsx")

BOYCEvaluation <- data.frame()
  
# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Get the current species name
  
  # Read model projection data for the species
  proj_names <- readRDS(paste0("models/", Sp, "/run_data.RDS"))

  ############# 1. K-fold cross validation

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

  ############# 2. Model formatting

  # 2.1 Model tunning
  #### A. XGBOOST
  
  # After tuning it on a separate Rscript (XGBoost_param)

  # Find the row index for the species
  row_index_XGBOOST <- which(XGBOOST_param$species == Sp)

  user.XGBOOST <- list('_allData_allRun' = list(
    nrounds = XGBOOST_param[row_index_XGBOOST, "nrounds"],
    max.depth = XGBOOST_param[row_index_XGBOOST, "max_depth"],
    eta = XGBOOST_param[row_index_XGBOOST, "eta"]))

  #### B. Random forest
  # After tuning it on a separate Rscript (XGBoost_param)
  
  # Find the row index for the species
  row_index_RF <- which(RF_param$species == Sp)

  user.RF <- list('_allData_allRun' = list(
    ntrees = RF_param[row_index_RF, "ntrees"],
    mtry = RF_param[row_index_RF, "mtry"]))

  # 2.2 Tuned model combination
  user.val <- list(
     XGBOOST.binary.xgboost.xgboost = user.XGBOOST,
     RF.binary.randomForest.randomForest = user.RF
     )
  
  # 2.3 Model parametrisation for running
  
  ### Define modeling options
  model_parameters <- bm_ModelingOptions(
    data.type = "binary", # Binary classification (presence/absence)
    models = c(
      'RF', 
      'XGBOOST',
      'MAXNET',
      'GAM',
      'GLM'),  # Full list of models
    strategy = "user.defined",
    user.val = user.val,
    bm.format = proj_names, # Input model data
    calib.lines = table_cv # Cross-validation table
  )

  
  ############# 3. Model running

  # # # Run modeling with specified parameters
  myBiomodModelOut <- BIOMOD_Modeling(
    proj_names,
    modeling.id = "Modeling",
    models = c('RF', 'XGBOOST','MAXNET',
               'GAM',
               'GLM'
               ),
    OPT.strategy = "user.defined",
    OPT.user = model_parameters,
    CV.strategy = 'user.defined',
    CV.user.table = table_cv,
    CV.do.full.models = FALSE,
    var.import = 10, # Number of variable importance metrics to calculate
    metric.eval = c('BOYCE'), # Evaluation metric (Boyce index)
    do.progress = TRUE)

  # Save the modeling output
  saveRDS(myBiomodModelOut, file = paste0("models/", Sp, "/output_models.rds"))

  
  ## 3.1 Response curves
  
  ## Get the evaluation results for each model
  evals <- get_evaluations(myBiomodModelOut)
  

  # Save model evaluation for all species
  Eval_Model_df <- data.frame(
    Species = Sp,
    Pseudo_absences = evals$PA,
    Run = evals$run,
    Model = evals$algo,
    Evaluation = evals$validation
  )
  BOYCEvaluation <- rbind(BOYCEvaluation, Eval_Model_df)
  
  
  ## Filter models with a Boyce index greater than 0.7
  Choosen_Model <- evals  %>% 
    filter(validation > 0.7)
  

  ## Extract the names of selected models
  selected_models <- Choosen_Model$full.name # Get the full names of selected models
   
  # Save the selected models
  saveRDS(selected_models, file = paste0("models/", Sp, "/selected_models.rds"))

  
  # Plot response curves for the selected models
  resp <- bm_PlotResponseCurves(bm.out = myBiomodModelOut,
                              models.chosen = selected_models,
                              fixed.var = "mean", # Fix the variable at the mean value
                              data_species = proj_names@data.species,
                              do.plot = FALSE, # Do not plot here (only generate the response table)
                              do.progress = FALSE)$tab

  colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")

  # Extract the model type from the 'Model' column
  resp <- resp %>%
    mutate(Model_Type = sub(".*_(GLM|GAM|RF|XGBOOST|MAXNET)$", "\\1", Model))

  # Define custom labels for each facet
  custom_labels <- c(
    "hurs_min" = "hurs_min",
    "bio5" = "bio5",
    "bio6" = "bio6",
    "globalCropland_2010CE" = "croplands",
    "npp" = "npp"
  )

  response <- ggplot(resp, aes(x = Var.value, y = Response)) +
    geom_line(alpha = 0.2, aes(group = Model, color = Model_Type)) +  # Use Model_Type for color
    stat_smooth(color = "black") + # Trend line
    facet_wrap(~Variable, scales = "free_x", labeller = as_labeller(custom_labels)) +
    theme_bw() +  # Black-and-white theme
    ylim(0, 1.1) +
    xlab("Variable value") +
    scale_color_manual(values = c(
      'RF' = '#26c031',
      'XGBOOST' = '#38599f',
      'MAXNET' = '#db1010',
      'GLM' = 'darkgrey',
      'GAM' = 'orange'
    )) +  # Custom color scale
    labs(color = "Model") +
    theme(
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      strip.text = element_text(size = 25, face = "bold")  # Adjust size and style here
    )


  # Create the output directory for figures if it doesn't exist
  save_dir <- paste0("output/Relation_Environnement")
  if(!dir.exists(save_dir)) {
  dir.create(save_dir)
  }

  # Save the plot
  ggsave(paste0(save_dir, "/",Sp,".jpeg",sep= ""),
         response,
       dpi = 2000,
       bg = NULL,
       width = 11,
       height = 8.5,
       units = "in")


  ## 3.2 Variable importance

  gg_varimp <- get_variables_importance(myBiomodModelOut) 


  colnames(gg_varimp) <- c("id",
                           "PA.Run",
                           "CV.Run", "Model",
                           "Variable",
                           "VI.run",
                           "Variable.importance")

  var_imp <- ggplot(gg_varimp, aes(y = Variable, x = Variable.importance)) +
    geom_boxplot() +
    theme_bw() +
    ggtitle(Sp) +
    scale_color_manual(values = c(
      'RF' = '#26c031',
      'XGBOOST' = '#38599f',
      'MAXNET' = '#db1010',
      'GLM' = 'darkgray',
      'GAM' = 'orange'
      )) +  # Custom color scale
    geom_jitter (alpha = 0.2, aes(col = Model))+
    theme(axis.text.x = element_text(size = 20),
                  axis.text.y = element_text(size = 20))

  # Create the output directory for figures if it doesn't exist
  save_dir <- paste0("output/Variable_Importance")
  if(!dir.exists(save_dir)) {
    dir.create(save_dir)
  }

  # Save the plot
  ggsave(paste0(save_dir, "/",Sp,".jpeg",sep=""),
         var_imp,
         dpi = 2000,
         bg = NULL,
         width = 11,
         height = 8.5,
         units = "in")
  
}

write.xlsx(BOYCEvaluation, file = "output/BOYCEvaluation.xlsx")


#### 4. BOYCE evaluation visualization

Vect_Sp <- rev(c("Hermetia illucens", "Tenebrio molitor", "Acheta domesticus", 
                  "Alphitobius diaperinus", "Musca domestica", 
                  "Gryllodes sigillatus", "Locusta migratoria",
                  "Gryllus assimilis"))

BOYCEvaluation$Species <- factor(BOYCEvaluation$Species, levels = Vect_Sp)
BOYCEvaluation <- BOYCEvaluation[order(BOYCEvaluation$Species), ]


boxplot_sp <- ggplot(BOYCEvaluation, aes(y = Species, x = Evaluation)) +
    geom_boxplot() + 
    theme_bw()  +
  geom_vline(aes(xintercept = 0.7), linetype = "dashed", color = "black") + # Add dashed line
    scale_color_manual(values = c(
      'RF' = '#26c031',
      'XGBOOST' = '#38599f',
      'MAXNET' = '#db1010',
      'GLM' = 'darkgray',
      'GAM' = 'orange'
      )) +  # Custom color scale
    geom_jitter (aes(col = Model))

# # Save the plot 
ggsave("figures/BoyceIndexselected.jpeg",
       boxplot_sp,
       dpi = 2000,
       bg = NULL,
       width = 11,
       height = 8.5,
       units = "in")
