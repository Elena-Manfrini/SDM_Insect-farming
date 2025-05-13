rm(list=ls())

# Libraries
library(terra)
library(raster)
library(ggplot2)
library(viridis)
library(openxlsx)
library(biomod2)
library(dplyr)

# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

# Save Highest risk categories for all species 
Projection_class_sp.all <- data.frame()

# Loop over each species to process occurrence data
for (i in 1:7) {
  Sp <- Vect_Sp[[i]] # Get the current species name.
  
  # Load pre-trained model outputs for the species
  myBiomodModelOut <- readRDS(paste0("models/", Sp, "/output_models.RDS"))
  selected_models <- readRDS(paste0("models/", Sp, "/selected_models.RDS"))
  
  # Remove XGBOOST & GAM models for projection
  selected_models_f <- selected_models[!grepl("XGBOOST", selected_models) & !grepl("GAM", selected_models)]
  

  ############# 1. Projection of individual models

  proj_runs <- BIOMOD_Projection(
    bm.mod = myBiomodModelOut, # Input individual model.
    proj.name = paste0(Sp, "_Projection"), # Name for the projection.
    new.env = Rastack, # Environmental data for projection.
    models.chosen = selected_models_f, # Only use selected models.
    build.clamping.mask = TRUE,
    output.format = ".tif")
  
  # Load suitability projection raster
  Projection_raw <- rast(paste0("models", "/", Sp, "/", gsub(" ", ".", Sp),
                                "/proj_", Sp, "_Projection",
                                "/proj_", Sp, "_Projection_", gsub(" ", ".", Sp), ".tif"))
  
  # Ensure no negative values in the raster
  Projection_raw <- app(Projection_raw, function(x) {
    x[x < 0] <- 0 # Replace negative values with 0.
    return(x)
  })

  ############# 2. Suitability raster

  ### 2.1 Projection Mean : Ensemble model
  
  Projection_ens <- mean(Projection_raw) # Calculate mean suitability across layers : Ensemble model
  names(Projection_ens) <- c("Suitability")
  
  # Save data frame
  # Create species-specific directories
  save_dir <- paste0("output/", "Suitability")
  if (!dir.exists(save_dir)) {
    dir.create(save_dir)
  }
  
  # Save raster of suitability
  writeRaster(Projection_ens, filename = paste0(save_dir, "/", Sp, "_ens_mod.tif"))
  
  # Convert raster to data frame for plotting
  raster_df <- as.data.frame(Projection_ens, xy = TRUE)
  
  
  ### 2.2 Projection Standard deviation 
  
  Projection_sd <- app(Projection_raw, sd) # Calculate standard deviation across layers.
  names(Projection_sd) <- c("Standard_deviation")
  
  # Save raster of standard deviation
  writeRaster(Projection_sd, filename = paste0("output/Suitability/", Sp, "sd_ens_mod.tif"))
  
  # Convert standard deviation raster to data frame for plotting
  raster_sd_df <- as.data.frame(Projection_sd, xy = TRUE)
 
  ### 2.3 Plots

  ## 2.3.1 Suitability map : Figure S5, Panel A
  
  # Species occurrences to get coordinates
  Fin_occ_var <- read.xlsx(paste0("data/filtered_occurences/Occ&Var_final", Sp, ".xlsx"))
  
  # Map
  plot_raster <- ggplot() +
    geom_tile(data = raster_df, aes(x = x, y = y, fill = Suitability)) +
    scale_fill_viridis(name = "Suitability", option = "D", direction = 1) +
    coord_equal() +  # Keep proportions
    theme_minimal() +
    labs(title = "Suitability Raster",
         subtitle = Sp
    ) +
    geom_point(data = Fin_occ_var,
               aes(x = x, y = y),
               colour = "brown",
               size = 0.2) +
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
    paste0(save_dir, "/Plot_Raster.jpeg", sep = ""),
    plot_raster,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  
  ## 2.2.1 Standard deviation map : Figure S5, Panel B
  
  # Plot standard deviation raster
  plot_raster_sd <- ggplot() +
    geom_tile(data = raster_sd_df, aes(x = x, y = y, fill = Standard_deviation)) +
    scale_fill_viridis(name = "Standard deviation", option = "inferno", direction = 1) +
    coord_equal() +  # Keep proportions
    theme_minimal() +
    labs(title = "Standard deviation Raster",
         subtitle = Sp,
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
    paste0("figures/", Sp, "/Plot_Raster_Standard-deviation.jpeg", sep = ""),
    plot_raster_sd,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
    

  
  ############# 3. Risk raster
  
  ### 2.1 Risk classification : Figure S6, Panel A
    
  ## Extract occurrences coordinates associated with suitability estimation
  Occu_Suit <- terra::extract(Projection_ens,
                              Fin_occ_var[, 1:2],
                              ID = FALSE)
  
  # Calculate occurrences quantiles (5% and 25%)
  quantiles <- quantile(Occu_Suit$Suitability, probs = c(0.05, 0.25))
  
  # Plot density histogram with quantile lines
  quantiles_plot <- ggplot(Occu_Suit, aes(x = Suitability)) +
    geom_density(fill = "#E0FFFF") +
    geom_vline(xintercept = quantiles[1], color = "black", linetype = "dashed", size = 1) +
    geom_vline(xintercept = quantiles[2], color = "black", linetype = "dashed", size = 1) +
    theme_minimal() +
    labs(
      title = paste0("Density Histogram of Suitability with Quantile Lines for ", Sp),
      x = "Suitability",
      y = "Density"
    ) +
    annotate("text", x = quantiles[1], y = 0.01, label = "5%", color = "black", vjust = -0.5, angle = 90) +
    annotate("text", x = quantiles[2], y = 0.01, label = "25%", color = "black", vjust = -0.5, angle = 90) +
    theme(
      axis.title.x = element_text(size = 20, face = "bold"),  # X axis title size and bold
      axis.title.y = element_text(size = 20, face = "bold"),  # Y axis title size and bold
      axis.text.x = element_text(size = 15),  # X axis tick labels size
      axis.text.y = element_text(size = 15)   # Y axis tick labels size
    )
  
  # Save the plot
  ggsave(
    paste0("figures/", Sp, "/DensityplotClass.jpeg", sep = ""),
    quantiles_plot,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  
  ### 2.2 Risk classification map visualization : Figure S6, Panel B
      
  # Classify raster values based on quantiles
  # Risk categories :
  # 0 -> 5% = Low or Unknown
  # 5 -> 25% = Intermediate
  # 25 -> 100% = High
  
  Projection_class <- Projection_ens
  
  Projection_class[Projection_class$Suitability < (quantiles[["5%"]])] <- 0
  Projection_class[Projection_class$Suitability >= quantiles[["5%"]] &
                     Projection_class$Suitability < quantiles[["25%"]]] <- 0.5
  Projection_class[Projection_class$Suitability >= quantiles[["25%"]]] <- 1
  
  # Save the Projection_class models
  writeRaster(Projection_class, filename = paste0("output/Suitability/", Sp, "_class_ens_mod.tif"))
  
  # Convert classified raster to data frame for plotting
  Projection_Class <- as.data.frame(Projection_class, xy = TRUE)
  
    
  # Plot of map risk categories
  plot_raster_Class <- ggplot() +
    geom_tile(data = Projection_Class, aes(x = x, y = y, fill = factor(Suitability))) +
    scale_fill_manual(
      name = "Suitability",
      values = c("0" = "#F5E8C2", "0.5" = "#F19134", "1" = "darkred"),  # Choose preferred colors
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
    paste0("figures/", Sp, "/Plot_RasterClass.jpeg", sep = ""),
    plot_raster_Class,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )

  ############# 3. Highest risk categories
  # Collect all High risk of species

  Projection_class_sp <- Projection_ens
  Projection_class_sp[Projection_class_sp$Suitability < quantiles[["25%"]]] <- 0
  Projection_class_sp[Projection_class_sp$Suitability >= quantiles[["25%"]]] <- 1
  
  # Convert classified raster to data frame for plotting
  Projection_class_sp <- as.data.frame(Projection_class_sp, xy = TRUE)
  colnames(Projection_class_sp) <- c("x", "y", "Hightsuit")
  
  # Add a column to track the species name
  Projection_class_sp$species <- Sp
  
  # Append to Projection_class_sp.all
  Projection_class_sp.all <- rbind(Projection_class_sp.all, Projection_class_sp)
}


### 3.1 Plot of cumulated High risk for all species : Figure 2

# Sum high risk of all species
Nbsp_Highsuit <- Projection_class_sp.all %>%
  group_by(x, y) %>%
  summarize(Nb_sp = sum(Hightsuit), .groups = "drop")

Nbsp_Highsuit$Nb_sp <- as.factor(Nbsp_Highsuit$Nb_sp)

# Map
  HighSuit_allsp <- ggplot() +
    geom_tile(data = Nbsp_Highsuit, aes(x = x, y = y, fill = Nb_sp)) +
    scale_fill_brewer(
      name = "Number of species with high suitability",
      palette = "YlOrRd",  # Use the YlOrRd palette
      direction = 1        # Colors go from low to high values
    ) +
    coord_equal() +  # Keep proportions
    theme_minimal() +
    labs(title = "") +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_blank(),
      axis.text = element_text(color = "gray50"),
      legend.position = "right",
      legend.key.height = unit(1, "cm")
    )

  ggsave(
    paste0("figures/Plot_RasterClass_allSP.jpeg",sep=""),
    HighSuit_allsp,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  