rm(list=ls())
library(terra)
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(biomod2)
library(dplyr)


# var.intervals <- intervals

###################### Define computeEnvCombinations function ######################
computeEnvCombinations <- function(env.stack,
                                   var.intervals,
                                   plot = TRUE,
                                   vars.to.plot = 1:5)
{
  # Convert the raster stack values to a data frame
  combinations <- as.data.frame(values(Rastack[[c(names(Rastack))]]))
  
  # Get coordinates for all cells in the raster stack and remove NA cells
  all.xy <- as.data.frame(xyFromCell(Rastack, 1:ncell(Rastack)))
  all.xy <- all.xy[complete.cases(combinations), ]
  
  combinations <- combinations[complete.cases(combinations), ]
  
  # All possible combination in the environment
  comb.cat <- sapply(colnames(combinations), function(x, combs, seqs.)
    cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
    combs = combinations, seqs. = var.intervals)      
  
  message("Total number of cells with environmental conditions in the geographical space: ", nrow(combinations),
          "\nNumber of duplicated conditions: ", length(which(duplicated(comb.cat))),
          "\nNumber of unique cells (environmental space): ", nrow(comb.cat[-which(duplicated(comb.cat)), ]))
  
  
  #################################################### Visualise duplicated environmental conditions
  
  comb.cat_df <- as.data.frame(comb.cat)
  
  cat_coords_sub_df <- comb.cat_df %>%
    mutate(All_val = apply(., 1, paste, collapse = "_"))  %>%
    group_by(All_val) %>%
    mutate(duplicated_comb = n() > 1) %>%
    ungroup()
  
  ### Bind coords with all var values
  cat_coords <- cbind(all.xy,cat_coords_sub_df[,c("duplicated_comb")])
  
  # # Species occurrences to get coordinates
  Fin_occ_var <- read.xlsx(paste0("data/filtered_occurences/Occ&Var_final_20Hermetia illucens.xlsx"))
  #
  
  ### Worldwide
  dup_world <- ggplot() +
    geom_tile(data = cat_coords, aes(x = x, y = y, fill = duplicated_comb)) +
    # scale_fill_viridis(name = "Environmental combination", option = "A", direction = 1) +
    coord_equal() +  # Keep proportions
    theme_minimal() +
    labs(title = "Duplicated environmental combinasion"
    ) +
    geom_point(data =  Fin_occ_var,
               aes(x = x, y = y),
               colour = "brown",
               size = 0.5) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14, face = "italic"),
      axis.title = element_blank(),
      axis.text = element_text(color = "gray50"),
      legend.position = "right",
      legend.key.height = unit(1, "cm"))
  
  
  # Save the plot
  ggsave(
    paste0("output/Dup_Env-cond_20.jpeg",sep=""),
    dup_world,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  
  ##### India 
  # Filter the dataframe based on x and y ranges : INDIA
  cat_coords_sub_df <- cat_coords %>%
    filter(x >= 68.7, x <= 97.25, y >= 8.4, y <= 37.6)
  
  Fin_occ_var_India <-  Fin_occ_var %>%
    filter(x >= 68.7, x <= 97.25, y >= 8.4, y <= 37.6)
  
  ## Check occurrences characteristics
  summary(merge(cat_coords_sub_df, Fin_occ_var_India, by = c("x", "y")))
  
  dup_india <- ggplot() +
    geom_tile(data = cat_coords_sub_df, aes(x = x, y = y, fill = duplicated_comb)) +
    # scale_fill_viridis(name = "Environmental combination", option = "A", direction = 1) +
    coord_equal() +  # Keep proportions
    theme_minimal() +
    labs(title = "Duplicated environmental combinasion"
    ) +
    geom_point(data =  Fin_occ_var_India,
               aes(x = x, y = y),
               colour = "black",
               size = 0.5) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14, face = "italic"),
      axis.title = element_blank(),
      axis.text = element_text(color = "gray50"),
      legend.position = "right",
      legend.key.height = unit(1, "cm"))
  
  # Save the plot
  ggsave(
    paste0("output/Dup_Env-cond_INDIA&Hermetiaocc_20.jpeg",sep=""),
    dup_india,
    dpi = 500,
    bg = NULL,
    width = 15,
    height = 8.5,
    units = "in"
  )
  
  ### Plot occurrences along environmental gradient (variable values)
  var.occ_2_long <- Fin_occ_var_India[, -c(1:3)] %>%
      tidyr::pivot_longer(cols = everything(), names_to = "Variable", values_to = "Values")

    # Create a list to store each ggplot
    plot_list <- list()
    # Define a list of colors for each variable
    colors <- c("bio5" = "brown",
                "hurs_min" = "#2BDBCA",
                "npp" = "green4",
                "bio6" = "black",
                # "Human_pop_2000" = "#8494FF",
                # "Human_footprint" = "blue",
                "globalCropland_2010CE" = "#E68613")

    # Loop over each environmental variable and create a boxplot for each
    for (var in unique(var.occ_2_long$Variable)) {
      plot_obj <- ggplot(var.occ_2_long[var.occ_2_long$Variable == var, ], aes(x = Variable, y = Values)) +
        geom_boxplot(fill = colors[var], color = "black", outlier.shape = 16) +
        labs(x = NULL, y = NULL, title = var) +
        theme_minimal() +
        theme(axis.text.x = element_blank(),  # Remove x-axis labels for individual plots
              axis.title.x = element_blank(),
              axis.text.y = element_text(size = 10),  # Adjust size of y-axis labels
              axis.title.y = element_text(size = 12)) +
        scale_color_manual(values = colors)   +    # Use colors as specified
        guides(color = "none")  # Remove the color legend
      # Add the plot to the list
      plot_list[[var]] <- plot_obj
    }
    #
    # Arrange all the plots in a grid
    all_plots <- grid.arrange(grobs = plot_list, ncol = 3)
    # # Save the grid to a PNG file
    ggsave(paste0("output/Filt_occurrences_plot/Variable_response_Hermetia_Inde.png"), plot = all_plots, width = 10, height = 8, dpi = 300)
    
    
  
  #####################################################################################################
  
  # Calculate midpoint values for each environmental variable interval
  possible.combs <- lapply(var.intervals, function(x)
    data.frame(interval = cut(x, x, right = FALSE, include.lowest = TRUE),
               mids = midpoints(cut(x, x, right = FALSE, include.lowest = TRUE))))
  
  # Assign midpoint values for each environmental variable in the spatial grid
  all.env.pixels <- sapply(colnames(comb.cat), function(x, int.to.replace, replacing.values)
  {
    replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
  }, int.to.replace = comb.cat, replacing.values = possible.combs)
  
  # Remove duplicate cells in environmental space
  duplicated.cells <- which(duplicated(all.env.pixels))
  if(length(duplicated.cells)){
    all.env.pixels <- all.env.pixels[-duplicated.cells, ]
    all.xy <- all.xy[-duplicated.cells, ]
  }
  
  ##### Check same envcondition worldwide :
  
  all.env.pixels_df <- data.frame(all.env.pixels)
  
  
  all.env.pixels_df <- all.env.pixels_df <- all.env.pixels_df %>%
    mutate(All_val = apply(., 1, paste, collapse = "_")) %>%  # Combine all values into one column
    group_by(All_val) %>%  
    mutate(group_id = cur_group_id()) %>%  # Assign unique number to each unique "All_val"
    ungroup()
  
  midpoint_coord_id <- cbind(all.xy,all.env.pixels_df)
  
  
  
  # Optional 3D plot of the selected environmental variables
  if(plot)
  {
    require(rgl)
    rgl::plot3d(all.env.pixels[, vars.to.plot])
  }
  # Return a list containing detailed intervals, unique conditions, and coordinates
  return(list(detailed.intervals = possible.combs,
              unique.conditions.in.env = all.env.pixels,
              coords = all.xy))
}

# Function to calculate midpoints of intervals
midpoints <- function(x) {
  lower <- as.numeric(gsub(",.*", "", gsub("\\(|\\[|\\)|\\]", "", x)))
  upper <- as.numeric(gsub(".*,", "", gsub("\\(|\\[|\\)|\\]", "", x)))
  midpoint <- lower + (upper - lower) / 2 ### ????? can't we just do the mean ??????
  return(midpoint)
}

#########################################################################################


############# 1. Environmental Space
# Load environmental raster stack
Rastack <- rast("data/final_baseline_20.tif")

# Get values from the raster stack, including NA for missing values (e.g., ocean areas)
values <- values(Rastack) 
combinations <- as.data.frame(values)
combinations <- combinations[complete.cases(combinations), ] # Remove rows with NA

# intervals <- list(
#   bio5 <- seq(floor(min(combinations[, 1])), ceiling(max(combinations[, 1])), length.out = 5),
#   hurs_min <- seq(floor(min(combinations[, 2])), ceiling(max(combinations[, 2])), by = 2),
#   # hurs_mean <- seq(floor(min(combinations[, 3])), ceiling(max(combinations[, 3])), by = 1),
#   # hurs_range <- seq(floor(min(combinations[, 4])), ceiling(max(combinations[, 4])), by = 1),
#   npp <- seq(floor(min(combinations[, 3])), ceiling(max(combinations[, 3])),by = 100)
#   # croplands <- seq(floor(min(combinations[, 6])), ceiling(max(combinations[, 6])), by = 0.01)
#   )

intervals <- list(
  bio5 = seq(min(combinations[, 1]), max(combinations[, 1]), length.out = 80),
  hurs_min = seq(min(combinations[, 2]), max(combinations[, 2]), length.out = 80),
  npp = seq(min(combinations[, 3]), max(combinations[, 3]), length.out = 80),
  globalCropland_2010CE = seq(min(combinations[, 4]), max(combinations[, 4]), length.out = 80),
  bio6 = seq(min(combinations[, 5]), max(combinations[, 5]), length.out = 80)
  # ,
  # Human_footprint = seq(min(combinations[, 5]), max(combinations[, 5]), length.out = 80),
  # Human_pop_2000 = seq(min(combinations[, 6]), max(combinations[, 6]), length.out = 80)
)
  
names(intervals) <- names(Rastack)
  
# Optionally save intervals for each step size as a separate file
saveRDS(intervals, file = paste0("data/intervals.rds"))
  
envir.space <- computeEnvCombinations(
      env.stack = Rastack,
      var.intervals = intervals,
      plot = TRUE,
      vars.to.plot = 1:5
    )
    
# Save environmental space for each step size
saveRDS(envir.space, file = paste0("data/Environmental_Space_20.rds"))
