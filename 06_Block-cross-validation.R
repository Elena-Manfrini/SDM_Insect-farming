library(openxlsx)
library(blockCV)
library(biomod2)
library(dplyr)

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$x

# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Current species name
  
  proj_names <- readRDS(paste0("models/", Sp, "/run_data.RDS"))
  proj_names
  ##############################################################################
  
  table_cv <- bm_CrossValidation(
    bm.format = proj_names,
    strategy = "kfold",
    k = 5,
    nb.rep = 1,
    balance = "both") ### aleatoire --> Changer en geographique ? Environmental ?
  
  calib_summary <-
    summary(proj_names, calib.lines =  table_cv) %>%
    filter(dataset == "calibration") ## separation jeu de données en 5
  
  proj_names@data.env.var
  
  iwp <- (10^6)^(1 - occurrences) ## parametre biomod
  
  RF_param_list <- NULL
  XGBOOST_param_list <- NULL
  MAXNET_param_list <- NULL
  
  for (cvrun in 1:nrow(calib_summary)) {
    prNum <- calib_summary$Presences[cvrun]
    bgNum <- calib_summary$Pseudo_Absences[cvrun]
    
    wt <- ifelse(occurrences == 1, 1, prNum / bgNum)

  # Ajustement des paramètres pour Random Forest
  if (prNum < bgNum) {
    RF_param_list[[paste0("_", calib_summary$PA[[cvrun]], "_", calib_summary$run[[cvrun]])]] <- list(
      ntree = 2000,  # Augmentation du nombre d'arbres
      mtry = floor(sqrt(ncol(Projection_Class))),  # Ajustement de mtry pour capturer plus de variables
      sampsize = c("0" = prNum, "1" = prNum),
      replace = TRUE
    )
  } else {
    RF_param_list[[paste0("_", calib_summary$PA[[cvrun]], "_", calib_summary$run[[cvrun]])]] <- list(
      ntree = 2000,  # Augmentation du nombre d'arbres
      mtry = floor(sqrt(ncol(Projection_Class))),  # Ajustement de mtry pour capturer plus de variables
      sampsize = c("0" = bgNum, "1" = bgNum),
      replace = TRUE
    )
  }
  # Ajustement des paramètres pour XGBOOST
  XGBOOST_param_list[[paste0("_", calib_summary$PA[[cvrun]], "_", calib_summary$run[[cvrun]])]] <- list(
    nrounds = 1000,  # Augmentation des itérations pour une meilleure convergence
    eta = 0.05,  # Augmentation du taux d'apprentissage pour accélérer la convergence
    max_depth = 7,  # Augmentation de la profondeur pour capturer plus d'interactions
    subsample = 0.9,  # Utilisation d'une plus grande fraction des données
    objective = "binary:logistic",
    gamma = 1,  # Introduction de la complexité pour éviter le surapprentissage
    colsample_bytree = 0.8,  # Augmentation de la fraction des caractéristiques à chaque arbre
    min_child_weight = 1,
    weight = wt,
    verbose = 0
  )
  
  # Ajustement des paramètres pour MAXNET
  MAXNET_param_list[[paste0("_", calib_summary$PA[[cvrun]], "_", calib_summary$run[[cvrun]])]] <- list(
    l1_regularizer = 0.1,  # Ajout d'une légère régularisation L1
    l2_regularizer = 0.1,  # Ajout d'une légère régularisation L2
    use_sgd = TRUE,
    set_heldout = 0,
    verbose = TRUE
  )
}

Param_Models <- list(
  RF.binary.randomForest.randomForest = RF_param_list,
  XGBOOST.binary.xgboost.xgboost = XGBOOST_param_list,
  MAXNET.binary.maxnet.maxnet = MAXNET_param_list
)



model_parameters <- bm_ModelingOptions(
  data.type = "binary",
  models = c('RF', 'XGBOOST','MAXNET'),
  strategy = "user.defined",
  user.base = "default",
  user.val = list(
    RF.binary.randomForest.randomForest = RF_param_list,
    XGBOOST.binary.xgboost.xgboost = XGBOOST_param_list,
    MAXNET.binary.maxnet.maxnet = MAXNET_param_list),
  bm.format = proj_names,
  calib.lines = table_cv) 


# Modélisation
myBiomodModelOut <- BIOMOD_Modeling(
  run_data,
  modeling.id = "1",
  models = c('RF', 'XGBOOST','MAXNET'),
  OPT.strategy = "user.defined",
  OPT.user = model_parameters,
  CV.strategy = 'user.defined',
  CV.user.table = table_cv,
  CV.do.full.models = FALSE,
  var.import = 10,
  metric.eval = c('BOYCE'),
  do.progress = TRUE
)

# Sauvegarder le résultat de la modélisation
saveRDS(myBiomodModelOut, file = paste0(save_dir, VecSp[i], "_model_output.rds"))

# Obtenir les évaluations des modèles
evals <- get_evaluations(myBiomodModelOut)

# Filtrer les modèles avec un indice de Boyce > 0.7
Choosen_Model <- evals %>% filter(validation > 0.7)

# Récupération des projections correspondantes
selected_models <- Choosen_Model$full.name # Nom des modèles sélectionnés


resp <- bm_PlotResponseCurves(bm.out = myBiomodModelOut,
                              models.chosen = selected_models,
                              # show.variables = cur_vars,
                              fixed.var = "mean",
                              data_species = occurrences$Observed,
                              do.plot = FALSE,
                              do.progress = FALSE)$tab

colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")

# setwd(str_c("C:/Users/ElenaPC/Documents/Quentin_M2/models",Epoque[j],"/",Method[l],sep=""))

p<-ggplot() +
  geom_line(data = resp[grep('RF', resp$Model), ],
            aes(x = Var.value, y = Response, group = Model), col = "darkgreen",alpha = 0.2) +
  geom_line(data = resp[grep('XGBOOST', resp$Model), ],
            aes(x = Var.value, y = Response, group = Model), col = 'darkblue', alpha = 0.2) +
  geom_line(data = resp[grep('MAXNET', resp$Model), ], 
            aes(x = Var.value, y = Response, group = Model), col = 'darkred', alpha = 0.2) +
  facet_wrap(~ Variable, scales = "free_x") +
  theme_bw() +
  ylim(0, 1.1) +
  xlab("Variable value")

ggsave(str_c("Relation_Environnement",VecSp[k],".jpeg",sep=""),
       p,
       # "jpeg",
       dpi = 2000,
       bg = NULL,
       width = 11,
       height = 8.5,
       units = "in"
)
}