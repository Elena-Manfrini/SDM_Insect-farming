rm(list = ls())

# Libraries
library(openxlsx)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)



############# 1. Species risk class preparation

# Upload Species name
Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$Vect_Sp

raster_list <- list() # Save class projection of all species

for (i in 1:7) {
  Sp <- Vect_Sp[[i]] # Get the current species name.
  
  # Read class projections
  Class_projection <- rast(paste0("output/Suitability/", Sp, "_class_ens_mod.tif"))
  
  raster_list[[Sp]] <- Class_projection
}

# Invasion risk of all species
Class <- rast(raster_list)
Class_df <- as.data.frame(Class, xy = TRUE)


############# 2. Insect farm data preparation

# Load insect farm dataset

Farms <- read_excel("data/Insect_farms.xlsx", 
                    sheet = "Final")

# Separate lat and long coordinates into 2 columns
Farms <- Farms %>%
  separate(`Coordonate Lat/Long`, into = c("y", "x"), sep = ",")

# As numeric coordinates
Farms$y <- as.numeric(Farms$y)
Farms$x <- as.numeric(Farms$x)

############# 3. Map visualization : Figure 5, Panel A

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


############# 4. Histogram : Figure 5, Panel B
# Percentage of farms in different risk categories

# Load environmental raster stack
Rastack <- rast("data/final_baseline.tif")

Order <- Farms[,c("x","y","Order")]

# Convert the data frame to a SpatVector
Order_vect <- vect(Order, geom = c("x", "y"), crs = crs(Rastack))

# Rasterize the farm locations, keeping the order information
Order_r <- rasterize(Order_vect, Rastack, field = "Order")

Order_r_df <- as.data.frame(Order_r, xy=T, na.rm = FALSE) # Convert to data frame with coordinates

SuitOrder <- merge(Order_r_df, Class_df, by = c('x','y'))

## 4.1 Diptera species formatting 

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

## 4.2 Coleoptera species formatting 

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

## 4.3 Orthoptera species formatting 

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

## 4.4 Histogram visualization

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
