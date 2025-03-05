rm(list=ls())
library(terra)
library(raster)
library(ggplot2)
library(viridis)
library(openxlsx)
library(biomod2)
library(dplyr)

# Load environmental raster stack
Rastack <- rast("data/final_baseline_20.tif")

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp


Projection_class_sp.all <- data.frame()


i <- 1

# Loop over each species to process occurrence data
for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]] # Get the current species name.
  
  # # Load pre-trained model outputs for the species
  myBiomodModelOut <- readRDS(paste0("models/", Sp, "/output_models_20_wtconv.RDS"))
  selected_models <- readRDS(paste0("models/", Sp, "/selected_models_20_wtconv.RDS"))

  ########################## 1. Projection of individual models

  proj_runs <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut, # Input individual model.
  proj.name = paste0(Sp,"_Projection"), # Name for the projection.
  new.env = Rastack, # Environmental data for projection.
  models.chosen = selected_models,# Only use selected models.
  build.clamping.mask = TRUE,
  output.format = ".tif")

  # Load suitability projection raster
  Projection_raw <- rast(paste0("models","/",Sp,"/",gsub(" ",".",Sp),
                                "/proj_",Sp,"_Projection",
                                "/proj_",Sp,"_Projection_",gsub(" ",".",Sp),".tif"))

  # Ensure no negative values in the raster
  Projection_raw <- app(Projection_raw, function(x) {
    x[x < 0] <- 0 # Replace negative values with 0.
    return(x)
  })
  
  ########################## 2. Suitability raster
  
   ### 2.1 Projection Mean : Ensemble model
  
  Projection_ens <- mean(Projection_raw) # Calculate mean suitability across layers : Ensemble model
  names(Projection_ens) <- c("Suitability")

  ## Save raster of suitability
  saveRDS(Projection_ens, file = paste0("output/Suitability/", Sp, "_ens_mod_20_wtconv.rds"))
  
  Projection_ens <- readRDS(paste0("output/Suitability/", Sp, "_ens_mod_20_wtconv.rds"))
 
   ###### Plot suitability raster
  
  # Convert raster to data frame for plotting
  raster_df <- as.data.frame(Projection_ens, xy = TRUE)
  
  
  #################################################################
  ### Test difference btw 0% and 5% convhull
  # ## Raster with 5%
  ens_20 <- readRDS(paste0("output/Suitability/", Sp, "_ens_mod_20.rds")) # Con 5%
  # Convert raster to data frame for plotting
  raster_df_20 <- as.data.frame(ens_20, xy = TRUE)
  colnames(raster_df_20) <- c('x', 'y', 'Suitability_CvHull')
  # 
  
  merge.Suit <- merge(raster_df,raster_df_20, by = c('x', 'y'))
  ## Hypothèse : davantage d'habitat favorables avec un convex hull qui elimine 2,5 % des outliers que 5%
  ## Différence entre les
  merge.Suit$Suitability_Difference <- merge.Suit$Suitability  - merge.Suit$Suitability_CvHull
  
  
  ####################### INDIA
  
  # Filter the dataframe based on x and y ranges : INDIA
  merge.Suit_India <- merge.Suit %>%
    filter(x >= 68.7, x <= 97.25, y >= 8.4, y <= 37.6)
  
  ### 4.1 Calculate quantiles (5% and 25%)
  
  Fin_occ_var <- read.xlsx(paste0("data/filtered_occurences/Occ&Var_final_20", Sp, ".xlsx"))
  
  ##### Without convhull
  Occu_Suit_wt <- cbind(Fin_occ_var[,1:2], 
                        terra::extract(ens_20, Fin_occ_var[,1:2], ID = FALSE))
  
  quantiles_wtconv <- quantile(Occu_Suit_wt$Suitability, probs = c(0.05, 0.25))
  nrow(merge.Suit_India[merge.Suit_India$Suitability < 484.446 , ]) * 100 / nrow(merge.Suit_India) 
  
  merge.Suit_India$Suitability[merge.Suit_India$Suitability < 484.446] <- 0
  merge.Suit_India$Suitability[merge.Suit_India$Suitability >= 484.446 & 
                                        merge.Suit_India$Suitability < 711.510] <- 0.5
  merge.Suit_India$Suitability[merge.Suit_India$Suitability > 711.510 ] <- 1
  
  ##### With convhull
  Occu_Suit <- terra::extract(Projection_ens,
                              Fin_occ_var[,1:2],
                              ID=F)
  quantiles_conv <- quantile(Occu_Suit$Suitability, probs = c(0.05, 0.25))
  nrow(merge.Suit_India[merge.Suit_India$Suitability_CvHull  < 568.006, ]) * 100 / nrow(merge.Suit_India) 
  
  merge.Suit_India$Suitability_CvHull[merge.Suit_India$Suitability_CvHull < 568.006 ] <- 0
  merge.Suit_India$Suitability_CvHull[merge.Suit_India$Suitability_CvHull >= 568.006 & 
                                        merge.Suit_India$Suitability_CvHull < 787.160] <- 0.5
  merge.Suit_India$Suitability_CvHull[merge.Suit_India$Suitability_CvHull >= 787.160 ] <- 1
  
  # Occurrences in India :
  Occu_Suit_India <- Occu_Suit_wt %>% # It does not matter if we choose with or without it's just coordinates which are interesting
    filter(x >= 68.7, x <= 97.25, y >= 8.4, y <= 37.6)
  
  
  ## PLot for visualisation
  ggplot() +
    # Highlight areas where Suitability == 0 in RED
    geom_tile(data = merge.Suit_India %>% filter(Suitability_CvHull == 0), 
              aes(x = x, y = y), fill = "red", alpha = 0.7) +
    # Highlight areas where Suitability == 1 in GREEN
    geom_tile(data = merge.Suit_India %>% filter(Suitability_CvHull == 0.5), 
              aes(x = x, y = y), fill = "orange", alpha = 0.7) +
    # Highlight areas where Suitability == 1 in GREEN
    geom_tile(data = merge.Suit_India %>% filter(Suitability_CvHull == 1), 
              aes(x = x, y = y), fill = "green", alpha = 0.7) +
    geom_point(data =  Occu_Suit_India,
               aes(x = x, y = y),
               colour = "black",
               size = 0.5)
  
  
  ######## Compare classes
  ## With conv hull
  Class20 <- readRDS(paste0("output/suitability/", Sp, "_class_ens_mod_20.rds"))
  Class20_df <- as.data.frame(Class20, xy = TRUE)
  colnames(Class20_df) <- c('x', 'y', 'Suitability_CvHull')
  # 
  Class_20wt <-readRDS(paste0("output/Suitability/", Sp, "_class_ens_mod_20_wtconv.rds"))
  Class_20wt_df <- as.data.frame(Class_20wt, xy = TRUE)
  colnames(Class_20wt_df) <- c('x', 'y', 'Suitability_20_Wt')
  # 
  merge.Suit_Class <- merge(Class20_df,Class_20wt_df, by = c('x', 'y'))
  # 
  # ## Hypothèse : davantage d'habitat favorables avec un convex hull qui elimine 2,5 % des outliers que 5%
  # ## Différence entre les
  merge.Suit_Class$SuitabilityClass_Difference <- merge.Suit_Class$Suitability_20_Wt - merge.Suit_Class$Suitability_CvHull
  merge.Suit_Class$SuitabilityClass_Difference <- as.factor(merge.Suit_Class$SuitabilityClass_Difference)
  
  table(merge.Suit_Class$SuitabilityClass_Difference)
  
  # Filter the dataframe based on x and y ranges : INDIA
  merge.Class_India <- merge.Suit_Class %>%
    filter(x >= 68.7, x <= 97.25, y >= 8.4, y <= 37.6)
  
  # 
  diff <- ggplot() +
    geom_tile(data = merge.Class_India, aes(x = x, y = y, fill = SuitabilityClass_Difference)) +
    scale_fill_manual(
      name = "Suitability changes",
      values = c("-0.5" = "darkred", "0" = "lightgrey", "0.5" = "darkgreen", "1" = "orange"),
      labels = c("Loose suitability", "No change", "Gain 0,5 suitability", "Gain 1 suitability")
    ) +
    coord_equal() +
    theme_minimal() +
    labs(
      title = "Suitability class differences between convexhull 0% - 5%",
      subtitle = Sp
    ) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14, face = "italic"),
      axis.title = element_blank(),
      axis.text = element_text(color = "gray50"),
      legend.position = "right",
      legend.key.height = unit(1, "cm")
    )
  # 
  # 
  # Save the plot
  ggsave(
    paste0("figures/", Sp,"/ClassDiff_20_0-5.jpeg",sep=""),
    diff,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  # 
  
  #################################################################
  ### Suitability map
  
  ## Save drataframe
  # Create species-specific directories
  save_dir <- paste0("output/", "Suitability")
  if (!dir.exists(save_dir)) {
    dir.create(save_dir)
  }

  # # # Save the dataframe of suitability
  # # saveRDS(raster_df, file = paste0("output/Suitability/", Sp, "_ens_mod_df_15_2,5.rds"))
  # # 
  # # Species occurrences to get coordinates
  Fin_occ_var <- read.xlsx(paste0("data/filtered_occurences/Occ&Var_final_20", Sp, ".xlsx"))
  # 
  # # Map
  plot_raster <- ggplot() +
    geom_tile(data = raster_df, aes(x = x, y = y, fill = Suitability)) +
    scale_fill_viridis(name = "Suitability", option = "D", direction = 1) +
    coord_equal() +  # Keep proportions
    theme_minimal() +
    labs(title = "Suitability Raster",
         subtitle = Sp
    ) +
    geom_point(data =  Fin_occ_var,
               aes(x = x, y = y),
               colour = "brown",
               size = 0.5) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14, face = "italic"),
      axis.title = element_blank(),
      axis.text = element_text(color = "gray50"),
      legend.position = "right",
      legend.key.height = unit(1, "cm"))


  # Create output directory for figures if it doesn't exist
  Save_dir_fig <- "figures"
  if (!dir.exists(Save_dir_fig)) {
    dir.create(Save_dir_fig)
  }

  # Create species-specific directories
  save_dir <- paste0("figures/", Sp)
  if (!dir.exists(save_dir)) {
    dir.create(save_dir)
  }

  # Save the plot
  ggsave(
    str_c("figures/", Sp,"/Plot_Raster_20_wtconv.jpeg",sep=""),
    plot_raster ,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  # 
  
  ########################## 3. Standard Deviation Raster
 
    Projection_sd <- app(Projection_raw, sd) # Calculate standard deviation across layers.
    names(Projection_sd) <- c("Standard_deviation")

    # Convert standard deviation raster to data frame for plotting
    raster_sd_df <- as.data.frame(Projection_sd, xy = TRUE)

    # Save the selected models
    saveRDS(raster_sd_df, file = paste0("output/suitability/", Sp, "_sd_ens_mod_20_wtconv.rds"))

    # Plot standard deviation raster
    plot_raster_sd <- ggplot() +
  geom_tile(data = raster_sd_df, aes(x = x, y = y, fill = Standard_deviation)) +
  scale_fill_viridis(name = "Standard deviation", option = "inferno", direction = 1) +
  coord_equal() +  # Pour garder les proportions
  theme_minimal() +  # Un thème propre
  labs(title = "Standard deviation Raster",
       subtitle = Sp,
       # caption = "Made with ggplot2"
       ) +
      theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "italic"),
    axis.title = element_blank(),
    axis.text = element_text(color = "gray50"),
    legend.position = "right",
    legend.key.height = unit(1, "cm"))

    # Save the standard deviation plot
    ggsave(
      str_c("figures/", Sp,"/Plot_Raster_Standard-deviation_20_wtconv.jpeg",sep=""),
      plot_raster_sd ,
      dpi = 500,
      bg = NULL,
      width = 15,
      height = 8.5,
      units = "in"
    )
  #   
  #   
  #   ##########################  4. Suitability classification
  # 
    Occu_Suit <- terra::extract(Projection_ens,
                          Fin_occ_var[,1:2],
                          ID=F)

    # boxplot(Occu_Suit)
    # boxplot(Occu_Suit[,"Suitability"]) # Boxplot of suitability values.

    ### 4.1 Calculate quantiles (5% and 25%)

    quantiles <- quantile(Occu_Suit$Suitability, probs = c(0.05, 0.25))
  #   
  #   # Plot density histogram with quantile lines
  #   ggplot(Occu_Suit, aes(x = Suitability)) +
  #     geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", alpha = 0.7, color = "black") +
  #     geom_density(color = "red", size = 1) +
  #     geom_vline(xintercept = quantiles[1], color = "black", linetype = "dashed", size = 1, show.legend = TRUE) +
  #     geom_vline(xintercept = quantiles[2], color = "black", linetype = "dashed", size = 1, show.legend = TRUE) +
  #     theme_minimal() +
  #     labs(
  #       title = "Density Histogram of Suitability with Quantile Lines",
  #       x = "Suitability",
  #       y = "Density"
  #     ) +
  #     annotate("text", x = quantiles[1], y = 0, label = "5%", color = "black", vjust = -0.5, angle = 90) +
  #     annotate("text", x = quantiles[2], y = 0, label = "25%", color = "black", vjust = -0.5, angle = 90)
  #    
    
    ### 4.2 Reclassify raster values based on quantiles
    
    ### 4.2.a 3 categories
    ## 0 -> 5% = Low
    ## 5 -> 25% = Intermediate
    ## 25 -> 100% = High
    
    Projection_class <-  Projection_ens

    Projection_class[Projection_class$Suitability < (quantiles[["5%"]])] <- 0
    Projection_class[Projection_class$Suitability >= quantiles[["5%"]] &
                       Projection_class$Suitability < quantiles[["25%"]]] <- 0.5
    Projection_class[Projection_class$Suitability >= quantiles[["25%"]]] <- 1

    # Save the selected models
    saveRDS(Projection_class, file = paste0("output/suitability/", Sp, "_class_ens_mod_20_wtconv.rds"))
    
     # Convert classified raster to data frame for plotting
    Projection_Class <- as.data.frame(Projection_class, xy=T)
    
    
  #   # Plot classified suitability raster
    plot_raster_Class <- ggplot() +
  geom_tile(data = Projection_Class, aes(x = x, y = y, fill = factor(Suitability))) +
  scale_fill_manual(
    name = "Suitability",
    values = c("0" = "#F5E8C2", "0.5" = "#F19134", "1" = "darkred"),  # Choisissez les couleurs que vous préférez
    labels = c("Low", "Intermediate", "High")
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Suitability capacity",
    subtitle = Sp
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "italic"),
    axis.title = element_blank(),
    axis.text = element_text(color = "gray50"),
    legend.position = "right",
    legend.key.height = unit(1, "cm"))

    ggsave(
      str_c("figures/", Sp,"/Plot_RasterClass_20_wtconv.jpeg",sep=""),
      plot_raster_Class,
      # "jpeg",
      dpi = 500,
      bg = NULL,
  width = 15,
  height = 8.5,
  units = "in"
  )
    

    ### 4.2.b Highest categories
    
#     Projection_class_sp <-  Projection_ens
#     Projection_class_sp[Projection_class_sp$Suitability < quantiles[["25%"]]] <- 0
#     Projection_class_sp[Projection_class_sp$Suitability >= quantiles[["25%"]]] <- 1
#     
#     # Convert classified raster to data frame for plotting
#     Projection_class_sp <- as.data.frame(Projection_class_sp, xy=T)
#     colnames(Projection_class_sp) <- c("x", "y", "Hightsuit")
#     
#     # Add a column to track the species name
#     Projection_class_sp$species <- Sp
#     
#     # Append to Projection_class_sp.all
#     Projection_class_sp.all <- rbind(Projection_class_sp.all, Projection_class_sp)
#     
 }
# 
# ############################## Plot combined High suitability data for all species
# ### **** 4.2.b Highest categories
# 
# Nbsp_Highsuit <- Projection_class_sp.all %>%
#   group_by(x, y) %>%
#   summarize(Nb_sp = sum(Hightsuit), .groups = "drop")
# 
#   
#   HighSuit_allsp <- ggplot() +
#   geom_tile(data = Nbsp_Highsuit, aes(x = x, y = y, fill = as.factor(Nb_sp))) +
#   scale_fill_viridis( discrete = TRUE, name = "Number of species with high suitability", option = "D", direction = 1) +  
#   coord_equal() +  # Keep proportions
#   theme_minimal() +  
#   labs(title = "") +
#   theme(
#     plot.title = element_text(size = 18, face = "bold"),
#     axis.title = element_blank(),
#     axis.text = element_text(color = "gray50"),
#     legend.position = "right",
#     legend.key.height = unit(1, "cm"))
# 
#   ggsave(
#     paste0("figures/Plot_RasterClass_allSP_20.jpeg",sep=""),
#     HighSuit_allsp,
#     # "jpeg",
#     dpi = 500,
#     bg = NULL, 
#     width = 15,
#     height = 8.5,
#     units = "in"
#   )
  