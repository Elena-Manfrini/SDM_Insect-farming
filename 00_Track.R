## Create all folder for clarity

library(openxlsx)
Vect_Sp <- c("Hermetia illucens", "Tenebrio molitor", "Acheta domesticus", 
             "Alphitobius diaperinus", "Musca domestica", 
             "Gryllodes sigillatus", "Locusta migratoria", "Gryllus assimilis")


# Data
dir.create("data")

dir.create("data/raw")
dir.create("data/raw/bioclim")
dir.create("data/raw/bioclim/baseline")
dir.create("data/raw/bioclim/land")
dir.create("data/raw/occurences")

dir.create("data/Filtered_occurences")

dir.create("data/ConvexHull")
for (i in Vect_Sp) {
  dir.create(paste0("data/ConvexHull/", i))
}

# model per species
dir.create("models")
dir.create("models/preparation")

for (i in 1:length(Vect_Sp)) {
  Sp <- Vect_Sp[[i]]
  dir.create(paste0("models/", Sp))
}

#For figures
dir.create("output")
dir.create("output/Filt_occurrences_plot")
