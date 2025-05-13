rm(list = ls())

# Libraries
library(readxl)
library(openxlsx)
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(httr) 
library(rnaturalearth)
library(ggplot2)
library(viridis)

# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

############# 1. IPBES Regions download

# Create IPBES folder if it does not exist
if(!dir.exists("data/IPBES_SubRegions")) {
  dir.create("data/IPBES_SubRegions")
}

### 1.1 Download IPBES file
# Uncomment the following lines to download and unzip the IPBES file
# adresse <- "https://zenodo.org/records/5719431/files/ipbes_regions_subregions_shape_1.1.zip?download=1"
# download.file(adresse, destfile = "data/IPBES_SubRegions/ipbes_subregions.zip", mode = "wb")
# unzip("data/IPBES_SubRegions/ipbes_subregions.zip", exdir = "data/IPBES_SubRegions")

### 1.2 Match IPBES and Raster format

# Load IPBES shapefile
IPBES_raw <- st_read("data/IPBES_SubRegions/IPBES_Regions_Subregions2.shp")

# Transform IPBES CRS to match the raster CRS
IPBES_CRS <- st_transform(IPBES_raw, crs = st_crs(Rastack[[1]]))

# Convert to spatial vector and crop to raster extent
IPBES_vect <- vect(IPBES_CRS)
IPBES_crop <- crop(IPBES_vect, Rastack[[1]])

# Define attribute columns to rasterize
attribute_columns <- c("ISO_3", "Area", "Sub_Region")

# Create a blank raster template
shape_raster <- rast(Rastack[[1]])

# Rasterize each attribute column and stack them
IPBES_list <- lapply(attribute_columns, function(col) {
  rasterize(IPBES_crop, shape_raster, field = col)
})

# Create a SpatRaster stack
IPBES_stack <- rast(IPBES_list)

# Save the final raster stack
writeRaster(IPBES_stack, filename = "data/IPBES.tif", overwrite = TRUE)

### 1.3 IPBES Regions & Subregions

IPBES_stack <- rast("data/IPBES.tif")

# Convert raster to data frame for IPBES SubRegions
IPBES_df <- as.data.frame(IPBES_stack, xy = TRUE)
IPBES_df <- IPBES_df[, -c(3)]

# Create a new column "Regions" based on "Sub_Region"
IPBESReg <- IPBES_df %>%
  mutate(Regions = case_when(
    grepl("america", Sub_Region, ignore.case = TRUE) | Sub_Region == "Caribbean" ~ "Americas",
    grepl("africa", Sub_Region, ignore.case = TRUE) ~ "Africa",
    grepl("asia", Sub_Region, ignore.case = TRUE) | Sub_Region == "Oceania" ~ "Asia and the Pacific",
    grepl("europe", Sub_Region, ignore.case = TRUE) | Sub_Region == "North East Asia" ~ "Europe and Central Asia",
    TRUE ~ NA_character_
  ))


############# 2. Figure preparation

IPBES_unique <- unique(IPBESReg[, c(3, 4, 5)])

# Alphabetic order
# Step 1: Sort by the first column (Area)
IPBES_unique <- IPBES_unique %>%
  arrange(Regions)

# Step 2: Sort by the second column (Sub_Region) within the fixed order of the first column
IPBES_unique <- IPBES_unique %>%
  arrange(Regions, Sub_Region)

# Step 3: Sort by the third column (Regions) within the fixed order of the first two columns
IPBES_unique <- IPBES_unique %>%
  arrange(Regions, Sub_Region, Area)

IPBES_unique <- IPBES_unique %>%
  arrange(3, 2, 1)

## Dataframes for countries
IPBESCountry_Perc <- IPBES_unique # Store percentage of the three different risk class per country
IPBESCountry_Freq <- IPBES_unique # Store the most frequent risk class per country

## Dataframe for subregions
IPBESubReg_all <- IPBES_unique # Store percentage of the three different risk class per subregion

## Dataframe for regions
IPBESReg_all <- IPBES_unique # Store percentage of the three different risk class per region

# Upload Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

Species_class <- data.frame() # Store the number of country in each risk category per species 

raster_list <- list() # Save class projection of all species

# Loop over each species to process occurrence data
for (i in 1:7) {
  Sp <- Vect_Sp[[i]] # Get the current species name.
  
  # Read class projections
  Class_projection <- rast(paste0("output/Suitability/", Sp, "_class_ens_mod.tif"))
  
  raster_list[[Sp]] <- Class_projection
  
  Class_projection <- as.data.frame(Class_projection, xy = TRUE)
  
  # Merge IPBES and Class projection by coordinates
  Class_IPBES <- merge(Class_projection, IPBESReg, by = c("x", "y"), all.y = TRUE)
  
  # Remove coordinates
  Class_IPBES <- Class_IPBES[,-c(1:2)]
  
  # Ensure Risk is treated as a character (not a factor)
  Class_IPBES$Suitability <- as.character(Class_IPBES$Suitability)
  
  
  ### 2.1 Data preparation for country
  ## 2.1.1 Percentage of risk categories per country
  
  # Calculate total count per area
  Perc_suit_Country <- Class_IPBES %>%
    group_by(Area, Suitability) %>%
    summarise(Count = n(), .groups = "drop") %>% ### Column which count the number of rows per suitablity category for each Area
    group_by(Area ) %>%
    mutate(Percentage = (Count / sum(Count)) * 100) %>% ### Percentage of each category per Area
    select(-Count) %>% # Remove the column
    pivot_wider(names_from = Suitability, values_from = Percentage, values_fill = 0) ## Suitability in column
  
  Perc_suit_Country <- Perc_suit_Country %>%
    rename_with(~ paste0(c("Low_risk_", "Intermediate_risk_", "High_risk_"), Sp), .cols = c("0", "0.5", "1"))
  
  # Dataframe of suitability class percentages per species per Subregions
  IPBESCountry_Perc <- merge(IPBESCountry_Perc, Perc_suit_Country, by = "Area")
  
  
  ## 2.1.2 Most frequent risk class per country
  # Attribute the most frequent risk class per country
  # Then, count the number of country in each category
  
  # Add a column indicating the most frequently occurring Suitability category.
  col_name <- paste0("Most_Frequent_Suitability_", Sp)
  
  Class_IPBES_Country_sp_freq <- Class_IPBES_Country_sp %>%
    group_by(Area) %>%
    mutate(!!col_name := factor(Suitability[which.max(tabulate(match(Suitability, unique(Suitability))))],
                                levels = levels(Suitability))) %>%
    ungroup()
  
  # Keep only one row per country
  unique_country <- unique(Class_IPBES_Country_sp_freq[,c(2,3)])
  
  # Dataframe of suitability the most frequent risk class per species per country
  IPBESCountry_Freq <- merge(IPBESCountry_Freq, unique_country, by = "Area")
  
  # Arrange for plotting 
  
  # Most frequent risk class per country
  colnames(unique_country) <- c("Area", "Most_Frequent_Suitability")
  
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
  
  # Number of county per suitability class per species
  Species_class <- rbind(Species_class,Nb_country_per_Suit_cat)
  
  
  
  ### 2.2 Data preparation for IPBES SubRegions
  ## Percentage of suitability classes per SubRegions
  
  ## Calculate total count per area
  Perc_suit_SubRegion <- Class_IPBES %>%
    group_by(Sub_Region, Suitability) %>%
    summarise(Count = n(), .groups = "drop") %>% ### Column which count the number of rows per suitablity category for each Area
    group_by(Sub_Region) %>%
    mutate(Percentage = (Count / sum(Count)) * 100) %>% ### Percentage of each category per Area
    select(-Count) %>% # Remove the column
    pivot_wider(names_from = Suitability, values_from = Percentage, values_fill = 0) ## Suitability in column
  
  Perc_suit_SubRegion <- Perc_suit_SubRegion %>%
    rename_with(~ paste0(c("Low_risk_", "Intermediate_risk_", "High_risk_"), Sp), .cols = c("0", "0.5", "1"))
  
  # Dataframe of suitability class percentages per species per Subregions
  IPBESubReg_all <- merge(IPBESubReg_all, Perc_suit_SubRegion, by = "Sub_Region")
  
  
  
  ### 2.3 Data preparation for IPBES Regions
  ## Percentage of suitability classes per Regions
  
  # Calculate total count per area
  Perc_suit_Region <- Class_IPBES %>%
    group_by(Regions, Suitability) %>%
    summarise(Count = n(), .groups = "drop") %>% ### Column which count the number of rows per suitablity category for each Area
    group_by(Regions) %>%
    mutate(Percentage = (Count / sum(Count)) * 100) %>% ### Percentage of each category per Area
    select(-Count) %>% # Remove the column
    pivot_wider(names_from = Suitability, values_from = Percentage, values_fill = 0) ## Suitability in column
  
  Perc_suit_Region <- Perc_suit_Region %>%
    rename_with(~ paste0(c("Low_risk_", "Intermediate_risk_", "High_risk_"), Sp), .cols = c("0", "0.5", "1"))
  
  # Add new column and assign name separately
  IPBESReg_all <- merge(IPBESReg_all, Perc_suit_Region, by = "Regions")
  
}

##### Save data frames
### Number of countries per suitability class for all species
Species_class <- as.data.frame(Species_class)

xlsx::write.xlsx(Species_class, "data/NbCountry-per-suitCat_All-sp.xlsx", row.names = F)

### Percentage of class category for all countries per species
xlsx::write.xlsx(IPBESCountry_Perc, "data/Perc-suitCat_percountry_All-sp.xlsx", row.names = F)

### Percentage of class category for all Subregions per species
xlsx::write.xlsx(IPBESubReg_all, "data/Perc-suitCat_perSubRegions_All-sp.xlsx", row.names = F)

### Percentage of class category for all Regions per species
xlsx::write.xlsx(IPBESReg_all, "data/Perc-suitCat_perRegions_All-sp.xlsx", row.names = F)


############# 3. IPBES figures

### 3.1 Figure 3 : Number of countries in each invasion risk category

Species_class <- read.xlsx("data/NbCountry-per-suitCat_All-sp.xlsx")

colnames(Species_class) <- sub(paste0("_risk"), "", colnames(Species_class))
Species_class <- Species_class[,c(1,4,3,2)]

reshaped_data <- Species_class %>%
  pivot_longer(
    cols = starts_with("High") | starts_with("Intermediate") |  starts_with("Low") ,
    names_to = "Risk_Level",
    values_to = "number"
  )

desired_orders <- rev(c("Hermetia illucens", "Tenebrio molitor", "Acheta domesticus", 
                        "Alphitobius diaperinus", "Musca domestica", 
                        "Gryllodes sigillatus", "Locusta migratoria"))

reshaped_data$Species <- factor(reshaped_data$Species, levels = desired_orders)
reshaped_data_ordered <- reshaped_data[order(reshaped_data$Species), ]


# Create the bar chart
histspnbcountry <- ggplot(reshaped_data, aes(x = Species, y = number, fill = Risk_Level)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  theme_minimal() +
  coord_flip() +
  scale_fill_manual(
    name = "Risk",
    values = c("High" = "darkred", "Intermediate" = "#F19134", "Low" = "#F5E8C2")
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Adjust size for x-axis text
    axis.text.y = element_text(size = 12),                       # Adjust size for y-axis text
    axis.title = element_text(size = 14),                         # Adjust size for axis titles
    legend.text = element_text(size = 12)                     # Adjust size for legend title
  ) +
  labs(
    x = "Species",
    y = "Number of country",
    fill = "Risk Level"
  )


ggsave(
  paste0("figures/NbCountryRisk.jpeg",sep=""),
  histspnbcountry,
  dpi = 500,
  bg = NULL,
  width = 15,
  height = 8.5,
  units = "in"
)

### 3.2 Figure 3 : Percentage of invasion risk of all IPBES sub-regions.

IPBESReg_all <- IPBESReg_all[,-c(1)]
IPBESReg_unique <- IPBESReg_all[!duplicated(IPBESReg_all$Sub_Region), ]

SubRegions <- read.xlsx("data/Perc-suitCat_perSubRegions_All-sp.xlsx")
SubRegions <- merge(SubRegions,IPBESReg_unique, by = "Sub_Region", all.x = TRUE)
# Step 1: Sort by the first column (Area)
SubRegions <- SubRegions %>%
  arrange(Regions)
# Step 2: Sort by the second column (Sub_Region) within the fixed order of the first column
SubRegions <- SubRegions %>%
  arrange(Regions, Sub_Region)


Vect_Sp <- c("Hermetia.illucens", "Tenebrio.molitor", "Acheta.domesticus", 
             "Alphitobius.diaperinus", "Musca.domestica", 
             "Gryllodes.sigillatus", "Locusta.migratoria")

Vect_Sp_file <- c("Hermetia illucens", "Tenebrio molitor", "Acheta domesticus", 
                  "Alphitobius diaperinus", "Musca domestica", 
                  "Gryllodes sigillatus", "Locusta migratoria")

for (i in  1:length(Vect_Sp)) {
  sp <- Vect_Sp[[i]]
  
  Sp <- Vect_Sp_file[[i]]
  
  # Create a subset of the dataframe
  subset_data <- SubRegions  %>%
    select(Sub_Region, contains(sp))
  
  colnames(subset_data) <- sub(paste0("_risk_", sp), "", colnames(subset_data))
  
  
  desired_orders <-  rev(c(
    "North America",
    "Caribbean",
    "Mesoamerica",
    "South America",
    "Eastern Europe",
    "Central and Western Europe",
    "North Africa",
    "Central Africa",
    "West Africa",
    "East Africa and adjacent islands",
    "Southern Africa",
    "Western Asia",
    "North-East Asia",
    "Central Asia",
    "South-East Asia",
    "South Asia",
    "Oceania"
  ))
  
  reshaped_data <- subset_data %>%
    pivot_longer(
      cols = starts_with("Low") | starts_with("Intermediate") | starts_with("High"),
      names_to = "Risk_Level",
      values_to = "Percentage"
    )
  
  
  reshaped_data$Sub_Region <- factor(reshaped_data$Sub_Region, levels = desired_orders)
  reshaped_data_ordered <- reshaped_data[order(reshaped_data$Sub_Region), ]
  
  reshaped_data_ordered$Sub_Region <- factor(reshaped_data_ordered$Sub_Region, levels = unique(reshaped_data_ordered$Sub_Region))
  
  
  # Create the stacked bar plot
  histSubReg <- ggplot(reshaped_data_ordered, aes(x = Sub_Region, y = Percentage, fill = Risk_Level)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    theme_minimal() +
    coord_flip() +
    scale_fill_manual(
      name = "Risk",
      values = c("Low" = "#F5E8C2", "Intermediate" = "#F19134", "High" = "darkred")
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Invasion Risk Levels of Hermetia illucens by SubRegions",
      x = "Country",
      y = "Percentage",
      fill = "Risk Level"
    )
  
  ggsave(
    paste0("figures/", Sp, "/SubRegionspercentage.jpeg",sep=""),
    histSubReg,
    dpi = 500,
    bg = NULL,
    width = 40,
    height = 40,
    units = "in"
  )
}


############# 3. Insect farms

### 3.1 Data preparation

## 3.1.1 Risk class

# Invasion risk of all species
Class <- rast(raster_list)
Class_df <- as.data.frame(Class, xy = TRUE)


## 3.1.1 Insect farm data preparation

# Load insect farm dataset
Farms <- read_excel("data/Insect_farms.xlsx", 
                    sheet = "Final")

# Separate lat and long coordinates into 2 columns
Farms <- Farms %>%
  separate(`Coordonate Lat/Long`, into = c("y", "x"), sep = ",")

# As numeric
Farms$y <- as.numeric(Farms$y)
Farms$x <- as.numeric(Farms$x)

### 3.2 Map visualization : Figure 5, Panel A

Category <- Farms[,c("x","y","Category")]

# Convert the data frame to a SpatVector
Farms_vect <- vect(Category, geom = c("x", "y"), crs = crs(Class))

# Rasterize the farm locations, keeping the order information
Farms_r <- rasterize(Farms_vect, Class, field = "Category")
Farms_r_df <- as.data.frame(Farms_r, xy=T, na.rm = FALSE) # Convert to data frame with coordinates

# Map preparation
raster_df <- Class_df[,c(1:3)]
raster_df$`Hermetia illucens` <- 1 # Unique value for a unique color (plot)

# Map
plot_farms <- ggplot() +
  geom_tile(data = raster_df, aes(x = x, y = y, fill = as.factor(`Hermetia illucens` ))) +
  scale_fill_manual(values = "lightgrey") +  # Set fill color to light grey
  coord_equal() +  # Keep proportions
  theme_minimal() +
  labs(title = "Insect farm location")  +
  geom_point(data = Farms_r_df,
             aes(x = x, y = y, color = Category, shape = Category),
             size = 3) +
  scale_color_manual(name = "Farm Category",
                     values = c("Orthoptera" = "#9A32CD", "Coleoptera" = "#E69F00", "Diptera" = "#00BFFF", "Multiple" = "darkred")) +
  scale_shape_manual(name = "Farm Order",
                     values = c("Orthoptera" = 16, "Coleoptera" = 17, "Diptera" = 18, "Multiple" = 15)) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_blank(),
    axis.text = element_text(color = "gray50"),
    legend.position = "none")


# Save the plot
ggsave(
  paste0("figures/Farmlocation.jpeg",sep=""),
  plot_farms,
  dpi = 500,
  bg = NULL,
  width = 15,
  height = 8.5,
  units = "in"
)


### 3.3 Histogram : Figure 5, Panel B
# Percentage of farms in different risk categories

Order <- Farms[,c("x","y","Order")]

# Convert the data frame to a SpatVector
Order_vect <- vect(Order, geom = c("x", "y"), crs = crs(Rastack))

# Rasterize the farm locations, keeping the order information
Order_r <- rasterize(Order_vect, Rastack, field = "Order")

Order_r_df <- as.data.frame(Order_r, xy=T, na.rm = FALSE) # Convert to data frame with coordinates

SuitOrder <- merge(Order_r_df, Class_df, by = c('x','y'))

## 3.3.1 Diptera species formatting 

Diptera_farms <- SuitOrder %>%
  filter(Order == "Diptera") %>%
  dplyr::select('Hermetia illucens', 'Musca domestica') %>%
  mutate(Diptera_Risk = pmax(`Hermetia illucens`, `Musca domestica`))

final_Dipt <- Diptera_farms %>%
  group_by(Diptera_Risk) %>%
  reframe(Order = "Diptera",
          FarmNumber = n(),
          Risk_Level = case_when(
            Diptera_Risk == 0.0 ~ "Low",
            Diptera_Risk == 0.5 ~ "Intermediate",
            Diptera_Risk == 1.0 ~ "High"
          )) %>%
  select(Order, FarmNumber, Risk_Level) %>%
  distinct()

## 3.3.2 Coleoptera species formatting 

Coleoptera_farms <- SuitOrder %>%
  filter(Order == "Coleoptera") %>%
  dplyr::select(`Tenebrio molitor`, `Alphitobius diaperinus`) %>%
  mutate(Coleoptera_Risk = pmax(`Tenebrio molitor`, `Alphitobius diaperinus`))

final_Col <- Coleoptera_farms %>%
  group_by(Coleoptera_Risk) %>%
  reframe(Order = "Coleoptera",
          FarmNumber = n(),
          Risk_Level = case_when(
            Coleoptera_Risk == 0.0 ~ "Low",
            Coleoptera_Risk == 0.5 ~ "Intermediate",
            Coleoptera_Risk == 1.0 ~ "High"
          )) %>%
  select(Order, FarmNumber, Risk_Level) %>%
  distinct()

## 3.3.3 Orthoptera species formatting 

Orthoptera_farms <- SuitOrder %>%
  filter(Order == "Orthoptera") %>%
  dplyr::select(`Gryllodes sigillatus`, `Acheta domesticus`, `Locusta migratoria`) %>%
  mutate(Orthoptera_Risk = pmax(`Gryllodes sigillatus`, `Acheta domesticus`, `Locusta migratoria`))

final_Orth <- Orthoptera_farms %>%
  group_by(Orthoptera_Risk) %>%
  reframe(Order = "Orthoptera",
          FarmNumber = n(),
          Risk_Level = case_when(
            Orthoptera_Risk == 0.0 ~ "Low",
            Orthoptera_Risk == 0.5 ~ "Intermediate",
            Orthoptera_Risk == 1.0 ~ "High"
          )) %>%
  select(Order, FarmNumber, Risk_Level) %>%
  distinct()

## 3.3.4 Histogram visualization

# Bind risk & farms for each order
Farmsrisk <- rbind(final_Dipt,final_Col)
Farmsrisk_final <- rbind(Farmsrisk,final_Orth)

# Percentage
Farmsrisk_final <- Farmsrisk_final %>%
  group_by(Order) %>%
  mutate(Percentage = (FarmNumber / sum(FarmNumber)) * 100)


histfarm <- ggplot(Farmsrisk_final, aes(x = Order, y = Percentage, fill = Risk_Level)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip()+
  scale_fill_manual(
    name = "Risk",
    values = c("Low" = "#F5E8C2", "Intermediate" = "#F19134", "High" = "darkred")
  ) +
  theme(axis.text.x = element_text(hjust = 1, , size = 20),  # Adjust size for x-axis text
        axis.text.y = element_text(size = 20),                       # Adjust size for y-axis text
        axis.title = element_text(size = 14),                         # Adjust size for axis titles
        legend.text = element_text(size = 12)) +
  labs(
    title = "Invasion Risk Levels by farms",
    x = "Species",
    y = "Number of farms",
    fill = "Risk Level"
  )

ggsave(
  paste0("figures/Farmrisk.jpeg",sep=""),
  histfarm,
  dpi = 500,
  bg = NULL,
  width = 15,
  height = 8.5,
  units = "in"
)
