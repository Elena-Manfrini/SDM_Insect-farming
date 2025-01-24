# Initial download and creation of environmental data
# Script to be run only once
rm(list=ls())
library(terra)
library(sf)
library(rnaturalearth)
library(rgbif)
library(openxlsx)


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
    options(timeout = 120)
  download.file(addresse,
                destfile = paste0("data/raw/bioclim/", bioclim, ".tif"),
                mode = "wb",
                quiet = TRUE)
}

#### 1.2 Environmental predictors cropland data from Cao et al., 2021
# We used cropland data from Cao et al., 2021 (https://doi.org/10.5194/essd-13-5403-2021)
# Download globalCropland_2010CE.tif and add it on "./data/donnees_brutes/bioclim" folder.

# Adresse of globalCropland_2010CE.tif
addresse_cropland <- "https://zenodo.org/records/5759237/files/globalCropland_2010CE.tif?download=1"
download.file(addresse_cropland,
              destfile = paste0("data/raw/bioclim/", "globalCropland_2010CE", ".tif"),
              mode = "wb",
              quiet = TRUE)

cropland <- c("globalCropland_2010CE") # Add cropland data
vars <- rbind(vars,cropland) # Add cropland variable to the previous variable list

#### 1.2 Human population density from Gao et al., 2020
# We used cropland data from Gao et al., 2020 (https://www.earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-pd-sspbsyr-1km-1.01)
# Download popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp1-geotiff.zip and add it on "./data/donnees_brutes/bioclim" folder.

human <- c("Human_pop_2000") # Add cropland data
vars <- rbind(vars,human) # Add cropland variable to the previous variable list

saveRDS(vars, "data/variable_names.RDS") # Save variable names as an RDS file
xlsx::write.xlsx(vars, "data/variable_names.xlsx", # Save variable names as an Excel file
                 row.names = FALSE)

############# 2. Species occurrences download (from GBIF)

Vect_Sp <- c("Hermetia illucens", "Tenebrio molitor", "Acheta domesticus", 
             "Alphitobius diaperinus", "Musca domestica", 
             "Gryllodes sigillatus", "Locusta migratoria", "Gryllus assimilis")

Vect_Sp.df <- as.data.frame(Vect_Sp)

xlsx::write.xlsx(Vect_Sp.df, "data/Species_names.xlsx",
                 row.names = FALSE)


log_file <- "data/raw/occurences/completion_log.txt"
output_dir <- "data/raw/occurences"


# Load completed species from the log file
completed_species <- if (file.exists(log_file)) {
  readLines(log_file)
} else {
  character()  # Initialize as empty if log file doesn't exist
}


for (i in 1:length(Vect_Sp)){
  Sp <- Vect_Sp[i] # Current species
  
  # Check if the species is already processed
  if (sp %in% completed_species) {
    message(sprintf("Skipping %s: already completed.", Sp))
    next
  }
  
  # Progress message
  message(sprintf("Processing %s (%d/%d)...", Sp, i, length(Vect_Sp)))
  
  
  # Download species occurrences from GBIF
  occurrence <- occ_search(scientificName = Sp,
                           basisOfRecord = c(
                             # "PRESERVED_SPECIMEN",
                           #                   "LIVING_SPECIMEN",
                           #                   "OBSERVATION",
                                             "HUMAN_OBSERVATION",
                                             "MACHINE_OBSERVATION"
                                             # "MATERIAL_SAMPLE",
                                             # "LITERATURE",
                                             # "MATERIAL_CITATION",
                                             # "OCCURRENCE"
                           ), ### remove "UNKNOWN","FOSSIL_SPECIMEN"
                           year = '1980,2024', 
                           limit = 300000,
                           fields= "minimal")
  
  Occu <- data.frame(
    x = c(occurrence$HUMAN_OBSERVATION$data$decimalLongitude, 
          occurrence$MACHINE_OBSERVATION$data$decimalLongitude
          # ,
          # occurrence$LIVING_SPECIMEN$data$decimalLongitude,
          # occurrence$OBSERVATION$data$decimalLongitude,
          # occurrence$MATERIAL_SAMPLE$data$decimalLongitude,
          # occurrence$LITERATURE$data$decimalLongitude,
          # occurrence$MATERIAL_CITATION$data$decimalLongitude,
          # occurrence$OCCURRENCE$data$decimalLongitude
          ),
    y = c(occurrence$HUMAN_OBSERVATION$data$decimalLatitude, 
          occurrence$MACHINE_OBSERVATION$data$decimalLatitude
          # ,
          # occurrence$LIVING_SPECIMEN$data$decimalLatitude,
          # occurrence$OBSERVATION$data$decimalLatitude,
          # occurrence$MATERIAL_SAMPLE$data$decimalLatitude,
          # occurrence$LITERATURE$data$decimalLatitude,
          # occurrence$MATERIAL_CITATION$data$decimalLatitude,
          # occurrence$OCCURRENCE$data$decimalLatitude
          )
  ) %>% na.omit()

  # Save occurrences to an Excel file
  xlsx::write.xlsx(Occu, paste0("data/raw/occurences/Occurences_Human&Machine_", Sp, ".xlsx"), row.names = FALSE)
  
  # Append the species to the completion log
  write(sp, file = log_file, append = TRUE)
  
  message(sprintf("Completed %s.", sp))
}
