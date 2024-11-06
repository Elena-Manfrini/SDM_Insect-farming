
library(terra)
library(sf)
library(Rarity) # plot des correlations multiples
library(virtualspecies)
library(sp)
library(geodata)
library(openxlsx)


############# 1. Rasterization of all variables

# Retrieve all the variable names
Vars <- read.xlsx("data/variable_names.xlsx")
Vect_Vars <- Vars$vars

# Create an empty list to store the rasters for each variable
raster_list <- list()
# Loop through each variable name and load the raster, converting to SpatRaster format if necessary
for (i in 1:length(Vect_Vars)) {
  Name_Var <- Vect_Vars[i] # Get the variable name
  # Load each raster file as a SpatRaster object
  raster_path <- paste0("data/raw/bioclim/", Name_Var, ".tif")
  raster <- rast(raster_path)
  # Add the raster to the list with its corresponding variable name
  raster_list[[Name_Var]] <- raster
}

# Ensure that all variables are on the same extent --> This is not the case here.
# Convert globalCropland_2010CE raster to same extent as CHELSA variables
raster_list[["globalCropland_2010CE"]] <- crop(raster_list[["globalCropland_2010CE"]], ext(raster_list[["CHELSA_npp"]]))

# stack rasters
Rastack <- rast(raster_list) # transform raster list into stack

r <- synchroniseNA(Rastack) ### from virtual species package

# Save the final raster stack with all variables
writeRaster(r, filename = "data/raw/bioclim/baseline.tif",overwrite = TRUE)

############# 2. Check colineratity

Rastack <- rast("data/raw/bioclim/baseline.tif")

#Create a PNG to store the output image of the collinearity tree
png("./output/collinearity_groups.png") # pour faire l'image avec l'arbre sur le disque

# Remove highly collinear variables by grouping them using a correlation threshold
groups <- removeCollinearity(raster::stack(Rastack), plot = T,
                             # This tests for collinearity and groups variables that are inter-collinear
                             multicollinearity.cutoff = 0.7, # Cutoff threshold
                             sample.points = TRUE,
                             nb.points = 20000,
                             method = "spearman") # Spearman correlation used due to possible non-normal variable distribution
dev.off()

############# 3. Construction of the final baseline environmental dataset

# Stack the definitive environmental variables into a single raster object
Rastack_fin <- Rastack[[c("CHELSA_bio5","CHELSA_bio7","CHELSA_hurs_min","CHELSA_hurs_range", "CHELSA_npp", "globalCropland_2010CE")]]

r_fin <- synchroniseNA(Rastack_fin) ### from virtual species package

## Save the final definitive raster
writeRaster(r_fin, filename = "./data/raw/bioclim/final_baseline.tif", overwrite = TRUE)

