rm(list = ls())
library(readxl)
library(openxlsx)
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(httr) # to download data off of Zenodo
library(rnaturalearth) # download ocean data from natural earth 
library(ggplot2)
library(viridis)

# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

# Create IPBES folder if it does not exist
if(!dir.exists("data/IPBES_SubRegions")) {
  dir.create("data/IPBES_SubRegions")
}


################### 1. IPBES Regions 

### 1.a Download IBPES file
# adresse <- "https://zenodo.org/records/5719431/files/ipbes_regions_subregions_shape_1.1.zip?download=1"
# download <- download.file(adresse, destfile = "data/IPBES_SubRegions/ipbes_subregions.zip", mode = "wb")
# IPBES_folder <- "data/IPBES_SubRegions/ipbes_subregions.zip"
# unzip(IPBES_folder, exdir = "data/IPBES_SubRegions")


### 1.b Match IPBES and Raster format

### IPBES
IPBES_raw <- sf::st_read("data/IPBES_SubRegions/IPBES_Regions_Subregions2.shp") # Shapefile

# Match CRS : WGS 84 (EPSG:4326)
IPBES_CRS <- st_transform(IPBES_raw, crs = st_crs(Rastack[[1]]))

# Crop and align the extent
# Convert to spat vector 
IPBES_vect <- vect(IPBES_CRS)

# Crop to raster extent
IPBES_crop <- crop(IPBES_vect, Rastack[[1]])

# Get column names (excluding geometry),
attribute_columns <- c("ISO_3", "Area","Sub_Region")

# Create a blank raster template with the same dimensions and resolution
shape_raster <- rast(Rastack[[1]])

# Rasterize each column and stack them
IPBES_list <- lapply(attribute_columns, function(col) {
  rasterize(IPBES_crop, shape_raster, field = col)
})

# Create a SpatRaster stack
IPBES_stack <- rast(IPBES_list)

# Save the final raster stack with all variables
writeRaster(IPBES_stack, filename = "data/IPBES.tif",overwrite = TRUE)




############## Figures ##############

IPBES_stack <- rast("data/IPBES.tif")

# Upload Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

### Convert raster to data frame
### IPBES Subregions
IPBESubReg <- as.data.frame(IPBES_stack[[3]], xy = TRUE)

### IPBES Countries
# Convert raster to data frame for plotting
IPBEScountry <- as.data.frame(IPBES_stack[[2]], xy = TRUE)

  
i <- 1
IPBESClass_all <- IPBESubReg
IPBEScountry_all <- IPBEScountry
Species_class <- data.frame()

# Loop over each species to process occurrence data
for (i in 1:7) {
  Sp <- Vect_Sp[[i]] # Get the current species name.
  
  
  
  #################### 2 Suitability  
  
  
  ## 2.a. Per Country
  
  ## Load raster of suitability
  Projection_ens <- readRDS(file = paste0("output/Suitability/", Sp, "_ens_mod_15.rds"))
  
  suit_df <- as.data.frame(Projection_ens, xy = TRUE)
  
  IPBESCountry_Suit <- merge(suit_df,IPBEScountry, by = c("x","y"))
  
  col_name <- paste0("SuitCountry_", Sp)
  
  mean_suit_by_country <- IPBESCountry_Suit %>%
    group_by(Area) %>%
    summarise(!!col_name := mean(Suitability, na.rm = TRUE))
  
  IPBES_SuitCountry_df <- IPBEScountry %>%
    left_join(mean_suit_by_country, by = "Area")
  
  
  plot_raster_Suit_Country <- ggplot() +
    geom_tile(data = IPBES_SuitCountry_df, aes(x = x, y = y, fill = !!sym(col_name))) +
    scale_fill_viridis(name = "Mean Suitability", option = "D", direction = 1) +  
    coord_equal() +  # Keep proportions
    theme_minimal() +  
    labs(title = "Suitability per countries",
         subtitle = Sp
    ) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14, face = "italic"),
      axis.title = element_blank(),
      axis.text = element_text(color = "gray50"),
      legend.position = "right",
      legend.key.height = unit(1, "cm"))
  
  # Save the plot
  ggsave(
    paste0("figures/", Sp,"/Suitability_per_country.jpeg",sep=""),
    plot_raster_Suit_Country ,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  
  
  
  ## 2.b IPBES regions
  
  IPBES_Sub_Suit <- merge(suit_df,IPBESubReg, by = c("x","y"))
  
  col_name <- paste0("Suitability_", Sp)
  
  mean_suit_by_Subregion <- IPBES_Sub_Suit %>%
    group_by(Sub_Region) %>%
    summarise(!!col_name := mean(Suitability, na.rm = TRUE))
  
  IPBES_Suit_df <- IPBESubReg %>%
    left_join(mean_suit_by_Subregion, by = "Sub_Region")
  
  
  plot_raster_Suit_Subregions <- ggplot() +
    geom_tile(data = IPBES_Suit_df, aes(x = x, y = y, fill = !!sym(col_name))) +
    scale_fill_viridis(name = "Mean Suitability", option = "D", direction = 1) +  
    coord_equal() +  # Keep proportions
    theme_minimal() +  
    labs(title = "Suitability Raster per IPBES Subregions",
         subtitle = Sp
    ) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14, face = "italic"),
      axis.title = element_blank(),
      axis.text = element_text(color = "gray50"),
      legend.position = "right",
      legend.key.height = unit(1, "cm"))
  
  # Save the plot
  ggsave(
    paste0("figures/", Sp,"/Suitability_per_Subregions.jpeg",sep=""),
    plot_raster_Suit_Subregions,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  
  
  
  #################### 3 Class projection
  
  
  # Read class projections
  Class_projection <- readRDS(file = paste0("output/suitability/", Sp, "_class_ens_mod_15.rds"))
  
  
  ##### 3.a Country
  
  ## Merge the two dataframes
  Class_IPBES_Country <- merge(Class_projection,IPBEScountry, by = c("x","y"))
  
  
  ##### 3.a.1. Number of country per category (low, intermediate, high) (Dataframe)
  
  # Keep only Countries and Suitability Classes
  Class_IPBES_Country_sp <- Class_IPBES_Country[,-c(1,2)]
  # As factor Suitability Classes
  Class_IPBES_Country_sp$Suitability <- as.factor(Class_IPBES_Country_sp$Suitability)
  
  # Add a column indicating the most frequently occurring Suitability category.
  Class_IPBES_Country_sp_test <- Class_IPBES_Country_sp %>%
    group_by(Area) %>%
    mutate(Most_Frequent_Suitability = factor(Suitability[which.max(tabulate(match(Suitability, unique(Suitability))))], 
                                              levels = levels(Suitability))) %>%
    ungroup()
  
  # Keep only one country row
  unique_country <- unique(Class_IPBES_Country_sp_test[,c(2,3)])
  
  # Number of country per suitability category
  Nb_country_per_Suit_cat <- unique_country %>%
    count(Most_Frequent_Suitability) %>%
    tidyr::pivot_wider(names_from = Most_Frequent_Suitability, values_from = n, values_fill = 0)
  
  Nb_country_per_Suit_cat <- Nb_country_per_Suit_cat %>%
    rename(
      "Low_risk" = `0`,
      "Intermediate_risk" = `0.5`,
      "High_risk" = `1`
    ) %>%
    mutate(Species = Sp) %>%
    select(Species, everything())  # Moves "Species" to the first column
  
  
  Species_class <- rbind (Species_class,Nb_country_per_Suit_cat)
  
  
  ##### 3.a.2. Map of High suitability percentage 
  
  ### Colname by species
  col_name <- paste0("PercHighSuit_", Sp)
  
  # Calculate percentage of High suitability per subregions
  percentage <- Class_IPBES_Country %>%
    group_by(Area) %>%
    summarise(!!col_name := sum(Suitability == 1, na.rm = TRUE) / n() * 100)
  
  # Merge percentage with initial IPBES dataframe
  IPBES_class_Country <- IPBEScountry %>%
    left_join(percentage, by = "Area")
  
  # Plot results
  plot_class_high_country <- ggplot() +
    geom_tile(data = IPBES_class_Country, aes(x = x, y = y, fill = !!sym(col_name))) +
    scale_fill_viridis(name = "Percentage of High Suitability", option = "C", direction = 1) +  
    coord_equal() +  # Keep proportions
    theme_minimal() +  
    labs(title = "High suitability class per IPBES country",
         subtitle = Sp
    ) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14, face = "italic"),
      axis.title = element_blank(),
      axis.text = element_text(color = "gray50"),
      legend.position = "right",
      legend.key.height = unit(1, "cm"))
  
  # Save the plot
  ggsave(
    paste0("figures/", Sp,"/HighSuitpercentage_per_country.jpeg",sep=""),
    plot_class_high_country ,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  
  
  ## 3.a.2 IPBES regions
  
  Class_IPBES_SuReg <- merge(Class_projection,IPBESubReg, by = c("x","y"))
  
  col_name <- paste0("PercHighSuit_", Sp)
  
  mean_suit_by_Subregion <- Class_IPBES_SuReg %>%
    group_by(Sub_Region) %>%
    summarise(!!col_name := mean(Suitability, na.rm = TRUE))
  
  IPBES_Suit_df <- IPBESubReg %>%
    left_join(mean_suit_by_Subregion, by = "Sub_Region")
  
  
  plot_raster_Class_Subregions <- ggplot() +
    geom_tile(data = IPBES_Suit_df, aes(x = x, y = y, fill = !!sym(col_name))) +
    scale_fill_viridis(name = "Percentage of High Suitability", option = "C", direction = 1) +  
    coord_equal() +  # Keep proportions
    theme_minimal() +  
    labs(title = "High suitability class per IPBES Subregions",
         subtitle = Sp
    ) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14, face = "italic"),
      axis.title = element_blank(),
      axis.text = element_text(color = "gray50"),
      legend.position = "right",
      legend.key.height = unit(1, "cm"))
  
  # Save the plot
  ggsave(
    paste0("figures/", Sp,"/HighSuitpercentage_per_Subregions.jpeg",sep=""),
    plot_raster_Class_Subregions,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  
  ## 3.a.3 Combine High suitability class of all species
  
  ## Bind with IPBESClass_F_all dataframe to combine percentage of all species
  final_IPBES <- IPBES_Suit_df[,-c(3)]
  IPBESClass_all <- merge(IPBESClass_all, final_IPBES, by = c("x","y"))
  
}

##### REF 3.a.1. Save data frame of Country number and Suitability category of all species

xlsx::write.xlsx(Species_class, paste0("data/NbCountry-per-suitCat_All-sp.xlsx"), row.names = F)




#################### 4. High suitability of all species

IPBESClass_all_mean <- IPBESClass_all %>%
  mutate(Mean_HighSuit = rowMeans(select(., starts_with("PercHighSuit"))))

IPBESClass_all_Fin <- IPBESClass_all_mean %>%
  mutate(Max_Suitability_Species = apply(select(., starts_with("PercHighSuit")), 1, function(x) names(x)[which.max(x)])) %>%
  mutate(Max_Suitability_Species = gsub("PercHighSuit_", "", Max_Suitability_Species))

## 4.1 Map of species with the highest suitability per subregions 

plot_Species<- ggplot() +
  geom_tile(data = IPBESClass_all_Fin, aes(x = x, y = y, fill = Max_Suitability_Species)) +
  scale_fill_manual( values = c("Hermetia illucens" = "#009688",
                     "Tenebrio molitor" =  "#673ab7",
                     "Acheta domesticus" = "#FFC107", 
                     "Alphitobius diaperinus" = "#FFA000",
                     "Musca domestica" = "#00796b",
                     "Gryllodes sigillatus" = "#ffecb3",
                     "Locusta migratoria" = "#ffeb3b",
                     "Gryllus assimilis" = "#9c27b0"),
    name = "Species") +
  # scale_fill_viridis(name = "Species with the highest suitability per Subregions", option = "C", direction = 1) +  
  coord_equal() +  # Keep proportions
  theme_minimal() +  
  labs(title = "Species with the highest suitability per Subregions"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "italic"),
    axis.title = element_blank(),
    axis.text = element_text(color = "gray50"),
    legend.position = "right",
    legend.key.height = unit(1, "cm"))

# Save the plot
ggsave(
  paste0("figures/SpeciesName_HighestSuitability.jpeg",sep=""),
  plot_Class_Subregions_allsp,
  dpi = 500,
  bg = NULL,
  width = 15,
  height = 8.5,
  units = "in"
)

## 4.2 Map of High suitability percentage mean of all species

plot_Class_Subregions_allsp <- ggplot() +
  geom_tile(data = IPBESClass_all_Fin, aes(x = x, y = y, fill = Mean_HighSuit)) +
  scale_fill_viridis(name = "Percentage of High Suitability", option = "C", direction = 1) +  
  coord_equal() +  # Keep proportions
  theme_minimal() +  
  labs(title = "Mean High suitability class per IPBES Subregions for all species"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "italic"),
    axis.title = element_blank(),
    axis.text = element_text(color = "gray50"),
    legend.position = "right",
    legend.key.height = unit(1, "cm"))

# Save the plot
ggsave(
  paste0("figures/HighSuitpercentage_per_Subregions_allsp.jpeg",sep=""),
  plot_Class_Subregions_allsp,
  dpi = 500,
  bg = NULL,
  width = 15,
  height = 8.5,
  units = "in"
)

## 4.3 Dataframe High suitability for all species

IPBES_SubR_Species <- IPBESClass_all_Fin[,-c(1,2)]
IPBES_SubR_Species.unique <- IPBES_SubR_Species[-which(duplicated(IPBES_SubR_Species)), ] 

xlsx::write.xlsx(IPBES_SubR_Species.unique, paste0("data/Per_per_SubRegion_All-sp.xlsx"), row.names = F)

  