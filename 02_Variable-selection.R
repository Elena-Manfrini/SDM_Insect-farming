library(terra)
library(sf)
library(Rarity) # plot des correlations multiples
library(virtualspecies)
library(sp)
library(geodata)
library(openxlsx)
library(rnaturalearth)
# library(rnaturalearthhires)
library(rnaturalearthdata)
library(ggplot2)
library(ggcorrplot) # for the correlation matrix

############# 1. Rasterization of all variables

# Retrieve all the variable names
Vars <- read.xlsx("data/variable_names.xlsx")
# Vars <- Vars[!Vars[, 1] %in% c("CHELSA_bio2", "CHELSA_bio7", "CHELSA_hurs_mean", "CHELSA_hurs_range", "globalCropland_2000CE"), ]
Vect_Vars <- Vars$vars

# Create an empty list to store the rasters for each variable
raster_list <- list()
# Loop through each variable name and load the raster, converting to SpatRaster format if necessary
for (i in 1:length(Vect_Vars)) {
  Name_Var <- Vect_Vars[i] # Get the variable name
  # Load each raster file as a SpatRaster object
  raster_path <- paste0("data/raw/bioclim/", Name_Var, ".tif")
  raster <- rast(raster_path)
  raster <- aggregate(raster, fact = 15, fun="mean") ## REMETTRE A 5: change from 1km^2 pixel to 5km^2
  # Add the raster to the list with its corresponding variable name
  raster_list[[Name_Var]] <- raster
}

# Upload Human pop density
human_pop <- rast("data/raw/bioclim/BaseYear_1km/baseYr_total_2000.tif")
human_pop <- aggregate(human_pop, fact = 5, fun="mean")
raster_list[[human_pop]] <- human_pop

# Ensure that all variables are on the same extent --> This is not the case here.
# Convert globalCropland_2010CE raster to same extent as CHELSA variables
raster_list[["globalCropland_2010CE"]] <- resample(raster_list[["globalCropland_2010CE"]],
                                                   raster_list[["CHELSA_npp"]])

# stack rasters
Rastack <- rast(raster_list) # transform raster list into stack

# To get values on continent and exclude thoose from oceans:
# Download land on natural earth: https://www.naturalearthdata.com/downloads/

# land_adresse <- "https://naciscdn.org/naturalearth/10m/physical/ne_10m_land.zip"
# download_land <- download.file(land_adresse, destfile = "data/raw/bioclim/land/land.zip")
# land_folder <- "data/raw/bioclim/land/land.zip"
# unzip(land_folder, exdir = "data/raw/bioclim/land")
# ne_download(scale = 10,
#             type = "land",
#             category = "physical",
#             destdir = "data/raw")
              
# Read the shapefile
land <- vect("data/raw/ne_10m_land.shp")

## Rastack with values only for land
Rastack <- mask(Rastack,land)

# Define the extent to crop: latitude between -60 and 75
extent_to_crop <- ext(-180.0001,  179.9999, -60, 75)

# Crop the raster stack
Rastack <- crop(Rastack, extent_to_crop)

# Save the final raster stack with all variables
writeRaster(Rastack, filename = "data/raw/bioclim/baseline.tif",overwrite = TRUE)

############# 2. Check collineratity
# Load Raster
Rastack <- rast("data/raw/bioclim/baseline.tif")

### 2.1 Correlation tree
#Create a PNG to store the output image of the collinearity tree
png("./output/collinearity_groups.png")

# Remove highly collinear variables by grouping them using a correlation threshold
groups <- removeCollinearity(Rastack, plot = T,
                             # This tests for collinearity and groups variables that are inter-collinear
                             multicollinearity.cutoff = 0.7, # Cutoff threshold
                             # sample.points = TRUE,
                             # nb.points = 50000,
                             method = "spearman") # Spearman correlation used due to possible non-normal variable distribution
dev.off()

##### Evaluation des variance inflation factor

vifs <- usdm::vif((Rastack)) # Toutes les variables

# Retenir une seule par groupe de variables intercorrélées
vifs <- usdm::vif(Rastack[[c("CHELSA_bio5", "CHELSA_hurs_min","CHELSA_npp")]])


### 2.2 Correlation matrix
# Rastack dataframe
Rastack_df <- values(Rastack)
Rastack_df <- as.data.frame(Rastack_df, xy = TRUE) # Keep only non-NA rows
Rastack_df <- Rastack_df[complete.cases(Rastack_df), ] ## remove rows with NA

#Calculate the correlation matrix
cor_matrix <- cor(Rastack_df, method = "pearson")

# plot the correlation matrix
#Create a PNG to store the output image of the collinearity tree
png("./output/correlation_matrix.png")
ggcorrplot(cor_matrix, method = "circle",  type = "upper")
dev.off()

############# 3. Construction of the final baseline environmental dataset

# Stack the definitive environmental variables into a single raster object
Rastack_fin <- Rastack[[c("CHELSA_bio5", "CHELSA_bio7", "CHELSA_hurs_min","CHELSA_hurs_range", "CHELSA_npp", "globalCropland_2010CE")]]

names(Rastack_fin) <- gsub("CHELSA_", 
                           "",
                           names(Rastack_fin))

## Save the final definitive raster
writeRaster(Rastack_fin, filename = "data/final_baseline.tif", overwrite = TRUE)
