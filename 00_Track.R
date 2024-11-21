## Create all folder for clarity

# Data
dir.create("data")
dir.create("data/raw")
dir.create("data/raw/bioclim")
dir.create("data/raw/bioclim/baseline")
dir.create("data/raw/bioclim/land")
dir.create("data/raw/occurences")
dir.create("data/Filtered_occurences")

# model per species
dir.create("models")
dir.create("models/preparation")

Species <- read.xlsx("data/Species_names.xlsx")
Vect_Sp <- Species$x

for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]]
  dir.create(paste0("models/", Sp))
}

#For figures
dir.create("output")
dir.create("output/Filt_occurrences_plot")
