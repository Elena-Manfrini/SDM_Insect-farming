# Initial download and creation of environmental data
# Script to be run only once

library(terra)
library(sf)
library(rnaturalearth)
library(rgbif)

############# 1. Environmental predictors download

#### 1.1 Environmental predictors from CHELSA
# 1.1.1 Create names for all the variables to be downloaded from CHELSA climate

# All variables are described in the pdf linked below
# https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf

vars <- c(paste0("bio", 1:11),  # Bioclimatic variables
          paste0("hurs_",  c("max", "mean", "min", "range")), # Surface relative humidity variables
          "npp" # Net primary productivity
)

vars <- data.frame(vars = vars) # Convert the variable names into a dataframe
vars$vars <- paste0("CHELSA_", vars$vars) # Prefix variables with "CHELSA_"

# 1.1.2  Download CHELSA data

for(i in 1:nrow(vars)) {
  bioclim <- vars[i,]
  bioclim_fin <- gsub("CHELSA_", "", bioclim)
    addresse <- paste0(
      "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_",
      bioclim_fin, "_1981-2010_", "V.2.1.tif")
  download.file(addresse,
                destfile = paste0("data/raw/bioclim/", bioclim, ".tif"),
                #method = "wget", 
                quiet = TRUE)
}

#### 1.2 Environmental predictors cropland data from Cao et al., 2021
# We used cropland data from Cao et al., 2021 (https://doi.org/10.5194/essd-13-5403-2021)
# Download globalCropland_2010CE.tif and add it on "./data/donnees_brutes/bioclim" folder.

# Adresse of globalCropland_2010CE.tif
addresse_cropland <- "https://zenodo.org/records/5105689/files/globalCropland_1000CE.tif?download=1"

download.file(addresse_cropland,
              destfile = paste0("data/raw/bioclim/", "globalCropland_2010CE", ".tif"),
              quiet = TRUE)

cropland <- c("globalCropland_2010CE") # Add cropland data
vars <- rbind(vars,cropland) # Add cropland variable to the previous variable list

saveRDS(vars, "data/variable_names.RDS") # Save variable names as an RDS file
xlsx::write.xlsx(vars, "data/variable_names.xlsx", # Save variable names as an Excel file
                  row.names = FALSE)

############# 3. Species occurrences download (from GBIF)

Species <- c("Hermetia illucens", "Tenebrio molitor", "Acheta domesticus", 
             "Alphitobius diaperinus", "Musca domestica", 
             "Gryllodes sigillatus", "Locusta migratoria", "Gryllus assimilis")

Species.df <- as.data.frame(Species)

xlsx::write.xlsx(Species, "data/Species_names.xlsx",
                 row.names = FALSE)

for (i in 1:length(Species)){
  sp <- Species[i] # Current species
  
  # Download species occurrences from GBIF
  occurrence <- occ_search(scientificName = sp, basisOfRecord = c("HUMAN_OBSERVATION","MACHINE_OBSERVATION"),
                           year = '1980,2024', limit = 300000, fields= "minimal")
  # Create a dataframe for species occurrences (latitude and longitude)
  Occu <- data.frame(
    x = c(occurrence$HUMAN_OBSERVATION$data$decimalLongitude, 
          occurrence$MACHINE_OBSERVATION$data$decimalLongitude),
    y = c(occurrence$HUMAN_OBSERVATION$data$decimalLatitude, 
          occurrence$MACHINE_OBSERVATION$data$decimalLatitude)
  ) %>% na.omit() # Remove rows with missing values (NA)
  #nrow(Occu) # Number of occurrences
  # Save occurrences to an Excel file
  xlsx::write.xlsx(Occu, paste0("data/raw/occurences/Occurences_", sp, ".xlsx"), row.names = F)
}
