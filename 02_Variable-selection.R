rm(list=ls())

# Libraries
library(terra)
library(sf)
library(Rarity)
library(virtualspecies)
library(sp)
library(geodata)
library(openxlsx)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(ggcorrplot) 

############# 1. Rasterization of all variables

# Retrieve all the variable names
Vars <- read.xlsx("data/variable_names.xlsx")
Vect_Vars <- Vars$vars

### 1.1 Bioclimatic variables

# Create an empty list to store the rasters for each variable
raster_list <- list()

# Loop through each variable name and load the raster, converting to SpatRaster format if necessary
for (i in 1:length(Vect_Vars)) {
  Name_Var <- Vect_Vars[i] # Get the variable name
  
  # Load each raster file as a SpatRaster object
  raster_path <- paste0("data/raw/bioclim/", Name_Var, ".tif")
  raster <- rast(raster_path)
  raster <- aggregate(raster, fact = 5, fun = "mean") # Change from 1km^2 pixel to 5km^2
  
  # Add the raster to the list with its corresponding variable name
  raster_list[[Name_Var]] <- raster
}

### 1.2 Human pop density
# Upload Human population density
Human_pop_2000 <- rast("data/raw/bioclim/Human_pop/baseYr_total_2000.tif")
Human_pop_2000 <- aggregate(Human_pop_2000, fact = 5, fun = "mean")
Human_pop_2000 <- resample(Human_pop_2000, raster_list[["CHELSA_npp"]])
raster_list[["Human_pop_2000"]] <- Human_pop_2000


### 1.3 Croplands
cropland <- rast("data/raw/bioclim/globalCropland_2010CE.tif")
cropland <- aggregate(cropland, fact = 5, fun = "mean")

# Ensure that all variables are on the same extent
# Convert globalCropland_2010CE raster to same extent as CHELSA variables
cropland <- resample(cropland, raster_list[["CHELSA_npp"]])
raster_list[["globalCropland_2010CE"]] <- cropland

### 1.4 Human footprint
# Upload Human footprint
Human_footprint <- rast("data/raw/bioclim/Human_footprint/wildareas-v3-2009-human-footprint.tif")
Human_footprint <- aggregate(Human_footprint, fact = 5, fun = "mean")

# Step 1: Reproject Human_footprint to match Human_pop_2000 CRS
Human_footprint_proj <- project(Human_footprint, Human_pop_2000)

Human_footprint <- resample(Human_footprint_proj, raster_list[["CHELSA_npp"]])
raster_list[["Human_footprint"]] <- Human_footprint_proj


# stack rasters
Rastack <- rast(raster_list) # transform raster list into stack


############# 2. Raster delimitation

# To get values on continent and exclude those from oceans:
# Download land on natural earth: https://www.naturalearthdata.com/downloads/

land_adresse <- "https://naciscdn.org/naturalearth/10m/physical/ne_10m_land.zip"
download_land <- download.file(land_adresse, destfile = "data/raw/bioclim/land/land.zip")
land_folder <- "data/raw/bioclim/land/land.zip"
unzip(land_folder, exdir = "data/raw/bioclim/land")
ne_download(scale = 10,
            type = "land",
            category = "physical",
            destdir = "data/raw")
              
# Read the shapefile
land <- vect("data/raw/ne_10m_land.shp")

## Rastack with values only for land
Rastack <- mask(Rastack, land)

# Define the extent to crop: latitude between -60 and 75
extent_to_crop <- ext(-180.0001, 179.9999, -60, 75)

# Crop the raster stack
Rastack <- crop(Rastack, extent_to_crop)

# Save the final raster stack with all variables
writeRaster(Rastack, filename = "data/raw/bioclim/baseline_20.tif", overwrite = TRUE)


############# 3. Check collinearity

Rastack <- rast("data/raw/bioclim/baseline.tif")

### 3.1 Correlation tree : Figure S6, Panel A

# Create a PNG to store the output image of the collinearity tree
png("./output/collinearity_groups.png")

# Remove highly collinear variables by grouping them using a correlation threshold
groups <- removeCollinearity(Rastack, plot = TRUE,
                             multicollinearity.cutoff = 0.6, # Cutoff threshold
                             method = "spearman") # Spearman correlation used due to possible non-normal variable distribution
dev.off()

### 3.2 Correlation matrix : Figure S6, Panel B

# Rastack dataframe
Rastack_df <- values(Rastack)
Rastack_df <- as.data.frame(Rastack_df, xy = TRUE) # Keep only non-NA rows
Rastack_df <- Rastack_df[complete.cases(Rastack_df), ] # Remove rows with NA

# Calculate the correlation matrix
cor_matrix <- cor(Rastack_df, method = "pearson")

# Plot the correlation matrix
# Create a PNG to store the output image of the collinearity tree
png("./output/correlation_matrix.png")
ggcorrplot(cor_matrix, method = "circle", type = "upper")
dev.off()

############# 4. Construction of the final baseline environmental dataset

# Stack the definitive environmental variables into a single raster object
Rastack_fin <- Rastack[[c("CHELSA_hurs_min", "CHELSA_bio5", "CHELSA_bio6", "globalCropland_2010CE", "CHELSA_npp")]]

names(Rastack_fin) <- gsub("CHELSA_", "", names(Rastack_fin))

## Save the final raster
writeRaster(Rastack_fin, filename = "data/final_baseline.tif", overwrite = TRUE)
