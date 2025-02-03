rm(list=ls())  # Clear workspace

# Load necessary libraries
library(terra)  # For raster and spatial operations
library(sf)  # For spatial features (optional for this script)
library(rgbif)  # For accessing GBIF data
library(openxlsx)  # For reading Excel files
library(taxize)  # For taxonomic information
library(raster)  # For raster operations
library(biomod2)  # For species distribution modeling

# Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

# Assign orders to species (Diptera, Coleoptera, Orthoptera)
Species$Order <- ifelse(Species$Vect_Sp %in% c("Hermetia illucens", "Musca domestica"), "Diptera", # ID 811
                        ifelse(Species$Vect_Sp %in% c("Tenebrio molitor", "Alphitobius diaperinus"), "Coleoptera", # ID 1470
                               ifelse(Species$Vect_Sp %in% c("Acheta domesticus", "Gryllodes sigillatus", "Locusta migratoria", "Gryllus assimilis"), "Orthoptera", # ID 1458  
                                      "Unknown")))
Vect_Order <- Species$Order


# Load environmental raster stack (representing environmental variables)
Rastack <- rast("data/final_baseline.tif")

# Change baseline Spat raster as raster
baseline_raster <- as(Rastack, "Raster")

# Transform raster as dataframe
Rastab <- as.data.frame(baseline_raster, xy=T, na.rm = FALSE)

i <- 3

for (i in 1:length(Vect_Order)) {
  Order <- Vect_Order[[i]] # Current order name
  
  ############# 1. Download order occurrences
  
  # 1.a  # Fetch taxonomic information using the taxon name
  order_info <- name_backbone(name = Order, rank = "order")
  # Extract the GBIF taxonKey (ID) for the order
  taxon_key_order <- order_info$usageKey
  
  # Subset species data by order
  sp.order <- subset(Species, Order == Vect_Order[[i]])
  
  # 1.b Download occurrence data for the order from GBIF
  occurrence.order <- occ_search(orderKey = taxon_key_order,
                         basisOfRecord = c(
                           "PRESERVED_SPECIMEN",
                           "LIVING_SPECIMEN",
                           "OBSERVATION",
                           "HUMAN_OBSERVATION",
                           "MACHINE_OBSERVATION",
                           "MATERIAL_SAMPLE",
                           "LITERATURE",
                           "MATERIAL_CITATION",
                           "OCCURRENCE"
                         ), ### remove "UNKNOWN","FOSSIL_SPECIMEN"
                         year = '1980,2024', 
                         limit = 5000,
                         fields= "minimal")
  
  # Format the occurrence data into a data frame
Occu.order <- data.frame(
  x = c(occurrence.order$HUMAN_OBSERVATION$data$decimalLongitude, 
        occurrence.order$MACHINE_OBSERVATION$data$decimalLongitude,
        occurrence.order$LIVING_SPECIMEN$data$decimalLongitude,
        occurrence.order$OBSERVATION$data$decimalLongitude,
        occurrence.order$MATERIAL_SAMPLE$data$decimalLongitude,
        occurrence.order$LITERATURE$data$decimalLongitude,
        occurrence.order$MATERIAL_CITATION$data$decimalLongitude,
        occurrence.order$OCCURRENCE$data$decimalLongitude
  ),
  y = c(occurrence.order$HUMAN_OBSERVATION$data$decimalLatitude, 
        occurrence.order$MACHINE_OBSERVATION$data$decimalLatitude,
        occurrence.order$LIVING_SPECIMEN$data$decimalLatitude,
        occurrence.order$OBSERVATION$data$decimalLatitude,
        occurrence.order$MATERIAL_SAMPLE$data$decimalLatitude,
        occurrence.order$LITERATURE$data$decimalLatitude,
        occurrence.order$MATERIAL_CITATION$data$decimalLatitude,
        occurrence.order$OCCURRENCE$data$decimalLatitude
        ),
  species = c(occurrence.order$HUMAN_OBSERVATION$data$scientificName, 
  occurrence.order$MACHINE_OBSERVATION$data$scientificName,
  occurrence.order$LIVING_SPECIMEN$data$scientificName,
  occurrence.order$OBSERVATION$data$scientificName,
  occurrence.order$MATERIAL_SAMPLE$data$scientificName,
  occurrence.order$LITERATURE$data$scientificName,
  occurrence.order$MATERIAL_CITATION$data$scientificName,
  occurrence.order$OCCURRENCE$data$scientificName)
  ) %>% na.omit()


############# 2. Preparing weigth map for background sampling

## 2.a Select unique value per species

# Remove space in scientificName column
Occu.order$species <- gsub(" ", "", Occu.order$species)
# Remove duplicated occurrences
Occu.order.unique <- Occu.order[-which(duplicated(Occu.order)), ] 

# Remove rows with missing coordinates or species names
if(length(is.na(which(is.na(Occu.order.unique$x))))) {
  Occu.order.unique <- Occu.order.unique[-which(is.na(Occu.order.unique$x)), ]
}


## 2.b Occurrences density preparation for background sampling
# Convert the occurrence dataframe into a spatial vector object 
Occu.order.unique_v <- vect(Occu.order.unique,
                 geom = c("x", "y"))
# Set the coordinate reference system (CRS) to WGS84 (EPSG:4326) (same as Rastack)
crs(Occu.order.unique_v) <- "EPSG:4326"

## 2.c Smooth Occurrences for defining the weight

density_order <- MASS::kde2d(crds(Occu.order.unique_v)[, 1],
                             crds(Occu.order.unique_v)[, 2],
                             n = dim(baseline_raster)[2:1],
                             lims = as.vector(ext(Rastack[[1]])))

rast_density_order = expand.grid(x = density_order$x, y = density_order$y, KEEP.OUT.ATTRS = FALSE)
rast_density_order$z = as.vector(density_order$z)
rast_density_order = rast(rast_density_order)

## 3. Shape density with land

# Read the shapefile
land <- vect("data/raw/ne_10m_land.shp")

## Keep density values only for land
rast_density_order <- mask(rast_density_order,land)

# Define the extent to crop: latitude between -60 and 75
extent_to_crop <- ext(-180.0001,  179.9999, -60, 75)

# Crop the raster stack
rast_density_order <- crop(rast_density_order, extent_to_crop)

# Ensures compatibility when resampling
rast_baseline_raster <- rast(baseline_raster)

# Resample the density raster to match the resolution and extent of the baseline raster
rast_density_order <- resample(rast_density_order,
                               rast_baseline_raster)

# Normalize the raster values by dividing each cell by the global maximum value
rast_density_order <- rast_density_order / global(rast_density_order, "max", na.rm = TRUE)[1, 1]

plot(rast_density_order)


## 4. Define background data
# Sample background points based on the density raster
background <- spatSample(rast_density_order,
                         method = "weights",
                         size = 2000,
                         replace = FALSE, # Pas de remise 
                         na.rm = TRUE, # Pas dans les données manquantes
                         xy = TRUE, # L'output inclut les coords XY
                         values = FALSE) # L'output exclut les variables

## 5. Prepare species data (occurrences + bias background)

j <- 1
for (j in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[j]] # Current species name
  Occu <- read.xlsx(paste0("data/raw/occurences/Occurences_", Sp, ".xlsx")) # Load occurrences for this species
  
  if (length(Occu$x) > 60){
    # Rasterize occurrences to align with baseline raster
    Occu_r <- rasterize(as.matrix(Occu[,1:2]), baseline_raster)
    Occu_r <- as.data.frame(Occu_r, xy=T, na.rm = FALSE) # Convert to data frame with coordinates
    
    Occu_r <- cbind(Occu_r, Rastab) # Add variable values
    Occu_r <- Occu_r[complete.cases(Occu_r), ] # Remove occurrences outside land and with no variable data
    Occu_r <- Occu_r[, c(1,2)]
    Occu_r[3] <- 1 # Assign presence = 1 for each occurrence
    colnames(Occu_r) <- c("x","y","Observed") # Rename columns
    
    # Add background points to occurrences
    P_points <- rbind.data.frame(Occu_r,
                                 data.frame(background, 
                                            Observed = 0))
    
    ### Visualisation data
    # Assign red for Observed = 1, and black for Observed = 0
    colors <- ifelse(P_points$Observed == 1, "red", "black")
    # Plot the raster
    plot(Rastack[[1]])
    # Plot the points with the conditional color
    points(P_points[, c("x", "y")], pch = 20, cex = 0.5, col = colors)
    
     
    # Create the output directory for models if it doesn't exist
    save_dir <- paste0("models/", Sp)
    if(!dir.exists(save_dir)) {
      dir.create(save_dir)
    }
    
    run_data <- BIOMOD_FormatingData(
      resp.name = Sp, # Species name
      resp.var = P_points[,3], # Présences + background
      expl.var = Rastack, # Explicative variables
      dir.name = save_dir, # Folder where we save models
      resp.xy = P_points[,c(1,2)], # xy coordinates of occurrences and background
      PA.strategy = NULL) # Pas de génération de points de background par biomod
    # Car on en a généré nous-mêmes
    
    saveRDS(run_data, file = paste0("models/", Sp, "/run_data_TGSP.RDS"))
  }
}
}