library(terra)
library(raster)
library(ggplot2)
library(viridis)
library(openxlsx)

# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

i <- 1
# Loop over each species to process occurrence data
for (i in 1:7) {
  Sp <- Vect_Sp[[i]] # Get the current species name.
  
  # Load pre-trained model outputs for the species
  myBiomodModelOut <- readRDS(paste0("models/", Sp, "/model_output_30itv.RDS"))
  selected_models <- readRDS(paste0("models/", Sp, "/selected_models_30itv.RDS"))
  
  ############# 1. Projection of Models
  
  ### 1.1 Individual models
  proj_runs <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut, # Input individual model.
  proj.name = paste0(Sp,"_Projection"), # Name for the projection.
  new.env = Rastack, # Environmental data for projection.
  models.chosen = selected_models,# Only use selected models.
  build.clamping.mask = TRUE,
  output.format = ".tif")
  
  ############# 2. Suitability raster
  
  # Load suitability projection raster
  Projection_raw <- rast(paste0("models","/",Sp,"/",gsub(" ",".",Sp),
                                "/proj_",Sp,"_Projection",
                                "/proj_",Sp,"_Projection_",gsub(" ",".",Sp),".tif"))
  
  # Ensure no negative values in the raster
  Projection_raw <- app(Projection_raw, function(x) {
    x[x < 0] <- 0 # Replace negative values with 0.
    return(x)
  })
  
  ### 2.1 Projection Mean : Ensemble model
  
  Projection_ens <- mean(Projection_raw) # Calculate mean suitability across layers : Ensemble model
  names(Projection_ens) <- c("Suitability")

  # Plotting suitability
  
  # Convert raster to data frame for plotting
  raster_df <- as.data.frame(Projection_ens, xy = TRUE)
  
  ## Save drataframe
  
  # Create species-specific directories
  save_dir <- paste0("output/", "Suitability")
  if (!dir.exists(save_dir)) {
    dir.create(save_dir)
  }
  
  # Save the selected models
  saveRDS(raster_df, file = paste0("output/Suitability/", Sp, "_ens_mod_30itv.rds"))
  
  ## Plot suitability raster
  
  # Occurrence of the species to get coordinates
  
  Fin_occ_var <- read.xlsx(paste0("data/filtered_occurences/Occ&Var_", Sp, ".xlsx"))
  
  # Map
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
    str_c("figures/", Sp,"/Plot_Raster_30itv.jpeg",sep=""),
    plot_raster ,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  
  ############# 3. Standard Deviation Raster
 
    Projection_sd <- app(Projection_raw, sd) # Calculate standard deviation across layers.
    names(Projection_sd) <- c("Standard_deviation")
    
    # Convert standard deviation raster to data frame for plotting
    raster_sd_df <- as.data.frame(Projection_sd, xy = TRUE)
    
    # Save the selected models
    saveRDS(raster_sd_df, file = paste0("output/suitability/", Sp, "_sd_ens_mod_30itv.rds"))
    
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
      str_c("figures/", Sp,"/Plot_Raster_Standard-deviation_30itv.jpeg",sep=""),
      plot_raster_sd ,
      dpi = 500,
      bg = NULL, 
      width = 15,
      height = 8.5,
      units = "in"
    )
    
    ############# 4. Suitability classification
    
    # Load occurrence data for the species
    Occu <- read.xlsx(paste0("data/filtered_occurences/Occ&Var_", Sp, ".xlsx"))
    
    # Extract suitability values at occurrence points
    Occu_Suit <- terra::extract(Projection_ens,
                          Occu[,1:2],
                          ID=F)
    
    boxplot(Occu_Suit) # Boxplot of suitability values.
  
    ### 4.1 Calculate quantiles (5% and 25%)
    
    quantiles <- quantile(Occu_Suit$Suitability, probs = c(0.05, 0.25))
    
    # Plot density histogram with quantile lines
    ggplot(Occu_Suit, aes(x = Suitability)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", alpha = 0.7, color = "black") +
      geom_density(color = "red", size = 1) +
      geom_vline(xintercept = quantiles[1], color = "black", linetype = "dashed", size = 1, show.legend = TRUE) +
      geom_vline(xintercept = quantiles[2], color = "black", linetype = "dashed", size = 1, show.legend = TRUE) +
      theme_minimal() +
      labs(
        title = "Density Histogram of Suitability with Quantile Lines",
        x = "Suitability",
        y = "Density"
      ) +
      annotate("text", x = quantiles[1], y = 0, label = "5%", color = "black", vjust = -0.5, angle = 90) +
      annotate("text", x = quantiles[2], y = 0, label = "25%", color = "black", vjust = -0.5, angle = 90)
    
    
    ### 4.2 Reclassify raster values based on quantiles
    ## 0 -> 5% = Low
    ## 5 -> 25% = Intermediate
    ## 25 -> 100% = High
    
     Projection_class <-  Projection_ens

     Projection_class[Projection_class$Suitability < (quantiles[["5%"]])] <- 0
     Projection_class[Projection_class$Suitability >= quantiles[["5%"]] & 
                        Projection_class$Suitability < quantiles[["25%"]]] <- 0.5
     Projection_class[Projection_class$Suitability >= quantiles[["25%"]]] <- 1

     # Convert classified raster to data frame for plotting
    Projection_Class <- as.data.frame(Projection_class, xy=T)
    write.table(Projection_Class,paste0("Projection_Class",Sp))
    
    # Plot classified suitability raster
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
      str_c("figures/", Sp,"/Plot_RasterClass_80itv.jpeg",sep=""),
      plot_raster_Class,
      # "jpeg",
      dpi = 500,
      bg = NULL, 
  width = 15,
  height = 8.5,
  units = "in"
  )
}











### High suitability per country 

# Convertir les données en un objet sf
Projection_sf <- st_as_sf(Projection_Class, coords = c("x", "y"), crs = 4326)

# Charger le shapefile des pays
world <- ne_countries(scale = "medium", returnclass = "sf")

# Rendre les géométries valides pour éviter les problèmes liés aux coordonnées non planaires
world_valid <- st_make_valid(world)

# Associer les pixels à un pays
pixels_countries <- st_join(Projection_sf, world_valid, join = st_intersects)

# Convertir en data frame pour éviter les problèmes de classe 'sf'
pixels_countries_df <- as.data.frame(pixels_countries)

# Calculer le nombre de pixels par classe de suitability et par pays
country_class_pixel_count <- pixels_countries_df %>%
  group_by(admin, Suitability) %>%
  summarize(pixel_count = n(), .groups = "drop")

# Calculer le nombre total de pixels par pays (toutes classes confondues)
total_pixels_per_country <- Projection_sf %>%
  st_join(world_valid, join = st_intersects) %>%
  as.data.frame() %>%
  group_by(admin) %>%
  summarize(total_pixels = n(), .groups = "drop")

# Joindre les données et calculer le pourcentage pour chaque classe
percentage_pixels <- left_join(country_class_pixel_count, total_pixels_per_country, by = "admin") %>%
  mutate(percentage = (pixel_count / total_pixels) * 100)

# Filtrer pour la classe de suitability "forte" (ou une autre classe spécifique)
percentage_pixels_strong <- percentage_pixels %>%
  filter(Suitability == 1)

# Convertir le shapefile en data frame pour éviter les problèmes avec 'left_join'
world_df <- as.data.frame(world)

# Joindre ce résultat au shapefile des pays
world_merged <- left_join(world_df, percentage_pixels_strong, by = c("admin" = "admin"))

# Convertir à nouveau en un objet sf pour la visualisation
world_merged_sf <- st_as_sf(world_merged)

# Créer la carte
carte_mondiale <- ggplot(data = world_merged_sf) +
  geom_sf(aes(fill = percentage), color = "black", size = 0.1) +
  scale_fill_viridis_c(
    name = "Pourcentage\nde surface par classe de suitability",
    na.value = "grey90",  # Pour les pays sans données
    option = "D",
    direction = 1
  ) +
  theme_minimal() +
  labs(
    title = "Pourcentage de surface globale par classe de suitability",
    subtitle = "pour chaque pays du monde",
    caption = "Source: Données raster et Natural Earth"
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "italic"),
    plot.caption = element_text(size = 10, hjust = 0),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_line(color = "gray70", size = 0.1)
  )




ggsave(
  
  str_c("Etablissement_par_pays",Sp,".jpeg",sep=""),
  
  carte_mondiale,
  
  # "jpeg",
  
  dpi = 500,
  
  bg = NULL,
  
  width = 15,
  
  height = 8.5,
  
  units = "in"
  
) 
