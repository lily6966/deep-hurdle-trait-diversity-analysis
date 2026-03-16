setwd("/Users/liyingnceas/GitHub/Deep-hurdle-population-model-exp/lightning_logs/historical")


#------------
library(fields)
require(rpart)
library(mgcv)
library(sf)
library(raster)
library(lubridate)
library(ranger)
library(scam)
library(PresenceAbsence)
library(verification)
library(fields)
library(gridExtra)
library(tidyverse)
library(sf)
library(raster)
library(lubridate)
library(ranger)
library(scam)
library(PresenceAbsence)
library(verification)
library(dggridR)
library(mlr)
library(dplyr)
library(ISOweek) 
library(rasterVis)
library(ggplot2)
library(terra)
library(readr)

library(tidyr)

library(readxl)
library(exactextractr)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection





#ttt

#--------------------------Presence prob mapping------------------------------------
#create directory to save result
tif_dir <- "occurance_plots_historical"
if (!dir.exists(tif_dir)) {
  dir.create(tif_dir)
}
pu <- vect("pu_multi_ccap_small.shp")
# get the buffered checklists for a given year and extract elevation values and calculate median and sd
regions <- pu %>%
  st_as_sf()


regions$group_id <- ceiling(seq_len(nrow(regions)) / 5)

regions_merged <- regions %>%
  group_by(group_id) %>%
  summarise()

plot(regions)
plot(regions_merged)



map_proj <- "ESRI:102003"
ne_land <- read_sf("gis-data-ccap.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry() %>%
  # project to the cdl projection
  st_transform(crs = 4326)

ne_country_lines <- read_sf("gis-data-ccap.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

# Get the bounding box of the trend_pu geometry
bbox <- st_bbox(regions)

# Center the trend_pu on the map
center_x <- (bbox$xmin + bbox$xmax) / 2
center_y <- (bbox$ymin + bbox$ymax) / 2


# Transform the bounding box to the same CRS as ne_land
bbox_buffered<-st_buffer(st_as_sfc(bbox), dist=50000)
bbox_transformed <- st_transform(bbox_buffered, crs = st_crs(bbox_buffered))

# Crop ne_land using the transformed bounding box
ne_land_cropped <- st_crop(ne_land, bbox_transformed)


predictions <- read_csv("presence_probability_results.csv") #1-12 pevelant species and 13-24 including rare species
species_num <- length(colnames(predictions))-2


for (j in 1:species_num){
  bird_species <- colnames(predictions)[j+2] # +2 to skip the first two columns (LON and LAT)
  # Create a SpatVector from coordinates
  points <- vect(predictions, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  # Define extent and resolution
  r_template <- rast(ext(points), resolution = 0.01, crs = crs(points))
  # Rasterize the column
  r_species <- rasterize(points, r_template, field = bird_species, fun = "mean")
  #clear trend dataframe
  trend_pu <- NULL
  #extract species-specific values from the raster
  pu_cell <- exact_extract(r_species, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(trend_mean = mean(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(regions, .)
  # bind to results
  
  trend_pu <- bind_rows(trend_pu, pu_cell)
  
  
  print(summary(trend_pu$trend_mean))
  # Clip trend_mean to a narrower range, e.g. -0.5 to 0.5 or -0.5 to 1.5
  #trend_pu$trend_scaled <- pmax(pmin(trend_pu$trend_mean, 1.5), -1)
  
  # Create the ggplot object for probability change
  abundance_change<-ggplot() +
    # Plot trend_pu with fill color based on trend_mean
    #geom_sf(data = trend_pu, aes(fill = trend_scaled), color = "#E2E2E2") +
    
    geom_sf(data = trend_pu, aes(fill = trend_mean), color = "#E2E2E2") +
    # Add base map
    geom_sf(data = ne_land_cropped, fill = "transparent", color = "#888888", lwd = 0.5) +
    # Add white layer for zero value in trend_mean
    
    
    #geom_sf(data = trend_pu[trend_pu$trend_mean == 0, ], fill = "white", color = "white") +
    
    scale_fill_distiller(
      palette = "YlGnBu",   # or "YlOrRd", "PuBu", etc.
      direction = 1,
      name = "Probability",
      limits = c(0, 1),
      na.value = "transparent"
    )+
    
    # scale_fill_gradient2(
    #   low = "#FFA500",    # orange
    #   mid = "white",
    #   high = "#228B22",   # forest green
    #   midpoint = 0,
    #   name = "Probability Change",
    #   limits = c(-1, 1.5),
    #   na.value = "transparent"
    # )
    
    # Gradient color scale for trend_mean
    #scale_fill_gradient(low = "#D7191C", high = "#3C7EA6", name = "Change in Relative Abundance", na.value = "transparent") +
    
    # Set plot title and caption
    labs(caption = paste0(bird_species)) +
    
    # Use minimal theme
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    # Center the trend_pu on the map
    coord_sf(crs = st_crs(ne_land_cropped), xlim = c(center_x - 2.5, center_x + 2.5), ylim = c(center_y-0.6, center_y+0.6)) 
  # Save the plot
  ggsave(paste0("occurance_plots_historical/", bird_species, "_occurance_probability_historical.pdf"), plot = abundance_change, width = 8, height = 6, units = "in")
}
   



#--------------------------Abundance mapping------------------------------------
#create directory to save result
plot_dir <- "abundance_plots_historical/"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)

}


predictions <- read_csv("multiplied_result.csv") #1-12 pevelant species and 13-24 including rare species
species_num <- length(colnames(predictions))-2
for (j in 1:species_num){
  bird_species <- colnames(predictions)[j+2] # +2 to skip the first two columns (LON and LAT)
  # Create a SpatVector from coordinates
  points <- vect(predictions, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  # Define extent and resolution
  r_template <- rast(ext(points), resolution = 0.01, crs = crs(points))
  # Rasterize the column
  r_species <- rasterize(points, r_template, field = bird_species, fun = "mean")
  
  #clear trend dataframe
  abundance_pu <- NULL
  #extract species-specific values from the raster
  pu_cell <- exact_extract(r_species, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(abundance_mean = mean(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(regions, .)
  # Skip this species if all values are NA or 0
  if (all(is.na(pu_cell$abundance_mean)) || all(pu_cell$abundance_mean %in% c(NA, 0))) {
    message(paste("Skipping", bird_species, "- only zeros or NAs"))
    next
  }
  # bind to results
  abundance_pu <- bind_rows(abundance_pu, pu_cell) #more useful for multiple years to stack
  
  
  
  # Create the ggplot object for probability change
  abundance_change<-ggplot() +
    # Plot trend_pu with fill color based on trend_mean
    geom_sf(data = abundance_pu, aes(fill = abundance_mean), color = "#E2E2E2") +
    # Add base map
    geom_sf(data = ne_land_cropped, fill = "transparent", color = "#888888", lwd = 0.5) +
    # Add white layer for zero value in trend_mean
    
    
    #geom_sf(data =range_pu[range_pu$range_mean == 0, ], fill = "white", color = "white") +
    
    scale_fill_distiller(
      palette = "YlOrRd",   # or "YlGnBu", "PuBu", etc.
      direction = 1,
      name = "Abundance",
      
      na.value = "transparent"
    )+
    
    # Set plot title and caption
    labs(caption = paste0(bird_species)) +
    
    # Use minimal theme
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    # Center the trend_pu on the map
    coord_sf(crs = st_crs(ne_land_cropped), xlim = c(center_x - 2.5, center_x + 2.5), ylim = c(center_y-0.6, center_y+0.6)) 
  
  ggsave(paste0(plot_dir, "/", bird_species, "_", "_abundance.pdf"), plot = abundance_change, width = 8, height = 6, units = "in")
  
  
}
# Install and load

library(TPD)
library(ks)
library(ggplot2)
library(gridExtra)

# === Step 1: Trait matrix ===
# Each row is a species, each column is a trait
count_predictions <-read_csv("multiplied_result.csv") #1-12 pevelant species and 13-24 including rare species
colSums(is.na(count_predictions))

occurance_predictions <-read_csv("presence_probability_results.csv")
count_predictions_inthigh <- read_csv("/Users/liyingnceas/GitHub/Deep-hurdle-population-model-exp/lightning_logs/inthigh/multiplied_result.csv") #1-12 pevelant species and 13-24 including rare species)
occurance_predictions_inthigh <- read_csv("/Users/liyingnceas/GitHub/Deep-hurdle-population-model-exp/lightning_logs/inthigh/presence_probability_results.csv")
# Convert species values to binary: 1 if > column mean, else 0
occurance_binary <- occurance_predictions %>%
  mutate(across(
    .cols = 3:ncol(.),
    .fns = ~ {
      col_mean <- mean(., na.rm = TRUE)
      if (col_mean < 0.01) {
        as.numeric(. > 0.05)
      } else {
        as.numeric(. > col_mean)
      }
    }
  ))

occurance_binary_inthigh <- occurance_predictions_inthigh %>%
  mutate(across(
    .cols = 3:ncol(.),
    .fns = ~ {
      col_mean <- mean(., na.rm = TRUE)
      if (col_mean < 0.01) {
        as.numeric(. > 0.05)
      } else {
        as.numeric(. > col_mean)
      }
    }
  ))
# 
count_filtered <- count_predictions
count_filtered[ , 3:ncol(count_filtered)] <- as.data.frame(
  ifelse(
    as.matrix(occurance_binary[ , 3:ncol(occurance_binary)]) == 1,
    as.matrix(count_predictions[ , 3:ncol(count_predictions)]),
    NA
  )
)

# Pivot to long format: one row per site-species pair filtered by species occurance
long_abundance <- count_filtered %>%
  pivot_longer(
    cols = -c(longitude, latitude),
    names_to = "species",
    values_to = "abundance"
  ) %>%
  filter(!is.na(abundance))

count_filted_inthigh <- count_predictions_inthigh
count_filted_inthigh[ , 3:ncol(count_filted_inthigh)] <- as.data.frame(
  ifelse(
    as.matrix(occurance_binary_inthigh[ , 3:ncol(occurance_binary_inthigh)]) == 1,
    as.matrix(count_predictions_inthigh[ , 3:ncol(count_predictions_inthigh)]),
    NA
  )
)
long_abundance_inthigh <- count_filtered_inthigh %>%
  pivot_longer(
    cols = -c(longitude, latitude),
    names_to = "species",
    values_to = "abundance"
  ) %>%
  filter(!is.na(abundance))

species_counts <- long_abundance %>%
  count(species, name = "n_obs")

species_counts_inthigh <- long_abundance_inthigh %>%
  count(species, name = "n_obs")

write_csv(species_counts, "species_obs_no.csv")


# Read trait file
traits_df <- read_excel("AVONET2_eBird.xlsx") %>%
  select(-Family2, -Order2, -Avibase.ID2, -Total.individuals, -Female, -Male, -Unknown, -Complete.measures, -Secondary1, -Traits.inferred,
         -Inference, -Reference.species) 

# Read taxonomy and map scientific to common names
taxonomy <- read_csv("ebird-taxonomy.csv")
sci_to_common <- setNames(taxonomy$common_name, taxonomy$scientific_name)
traits_df$Species2 <- sci_to_common[traits_df$Species2]
traits_df <- traits_df %>%
  drop_na(Species2) 
# Extract species names from count_predictions (exclude lat/lon)
species_names <- colnames(count_predictions)[-(1:2)]
# Filter out NA species
traits_filtered <- traits_df %>%
  filter(!is.na(Species2)) %>%
  # Keep only those species in the prediction data
  filter(Species2 %in% species_names)

# Step 1: Identify shared species
shared_species <- intersect(traits_filtered$Species2, species_names)

# Step 2: Subset `count_predictions` to keep only shared species
# Reorder to match `traits_filtered$Species2`
count_predictions_filtered <- count_predictions %>%
  select(all_of(c(traits_filtered$Species2)))
count_predictions_filtered_inthigh <- count_predictions_inthigh %>%
  select(all_of(c(traits_filtered$Species2)))

#remove non-trait columns
traits_filtered <- traits_filtered %>%
  select(-Mass.Refs.Other, -Mass.Source)
#merge trait means with observation counts
trait_obs_merged <- traits_filtered %>%
  rename(species = Species2) %>%
  inner_join(species_counts, by = "species")

trait_obs_merged_inthigh <- traits_filtered %>%
  rename(species = Species2) %>%
  inner_join(species_counts_inthigh, by = "species")



# Numeric trait columns to expand
# Identify numeric columns (excluding species column)
# Identify trait columns
numeric_cols <- traits_filtered %>%
  select(where(is.numeric)) %>%
  colnames()

categorical_cols <- traits_filtered %>%
  select(-all_of(numeric_cols), -Species2) %>%  # exclude numeric + species
  colnames()
# Expand traits per species
traits_expanded <- do.call(rbind, lapply(1:nrow(trait_obs_merged), function(i) {
  row <- trait_obs_merged[i, ]
  n_obs_row <- row$n_obs
  species <- row$species
  
  # Simulate numeric traits
  simulated_numeric <- lapply(numeric_cols, function(trait) {
    mu <- as.numeric(row[[trait]])
    if (is.na(mu) || mu == 0) {
      rep(NA, n_obs_row)
    } else {
      rnorm(n_obs_row, mean = mu, sd = abs(mu) * 0.05)
    }
  })
  
  # Repeat categorical traits
  repeated_categorical <- lapply(categorical_cols, function(trait) {
    rep(row[[trait]], n_obs_row)
  })
  
  # Combine everything
  simulated_df <- as.data.frame(c(simulated_numeric, repeated_categorical))
  colnames(simulated_df) <- c(numeric_cols, categorical_cols)
  
  simulated_df$species <- species
  return(simulated_df)
}))


# Expand traits per species
traits_expanded_inthigh <- do.call(rbind, lapply(1:nrow(trait_obs_merged_inthigh), function(i) {
  row <- trait_obs_merged_inthigh[i, ]
  n_obs_row <- row$n_obs
  species <- row$species
  
  # Simulate numeric traits
  simulated_numeric <- lapply(numeric_cols, function(trait) {
    mu <- as.numeric(row[[trait]])
    if (is.na(mu) || mu == 0) {
      rep(NA, n_obs_row)
    } else {
      rnorm(n_obs_row, mean = mu, sd = abs(mu) * 0.05)
    }
  })
  
  # Repeat categorical traits
  repeated_categorical <- lapply(categorical_cols, function(trait) {
    rep(row[[trait]], n_obs_row)
  })
  
  # Combine everything
  simulated_df <- as.data.frame(c(simulated_numeric, repeated_categorical))
  colnames(simulated_df) <- c(numeric_cols, categorical_cols)
  
  simulated_df$species <- species
  return(simulated_df)
}))




library(dplyr)

# 1. Arrange both dataframes by species name
long_abundance_ordered <- long_abundance %>%
  arrange(species)
long_abundance_ordered_inthigh <- long_abundance_inthigh %>%
  arrange(species)

traits_expanded_ordered <- traits_expanded %>%
  arrange(species)
traits_expanded_ordered_inthigh <- traits_expanded_inthigh %>%
  arrange(species)

# 2. Create a per-species row index
long_abundance_ordered <- long_abundance_ordered %>%
  group_by(species) %>%
  mutate(row_id = row_number()) %>%
  ungroup()
long_abundance_ordered_inthigh <- long_abundance_ordered_inthigh %>%
  group_by(species) %>%
  mutate(row_id = row_number()) %>%
  ungroup()


traits_expanded_ordered <- traits_expanded_ordered %>%
  group_by(species) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

traits_expanded_ordered_inthigh <- traits_expanded_ordered_inthigh %>%
  group_by(species) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

# 3. Join by species and row_id (ensures correct alignment)
traits_with_abund <- long_abundance_ordered %>%
  inner_join(traits_expanded_ordered, by = c("species", "row_id"))%>%
  #remove species with less than 4 observations
  group_by(species) %>%
  filter(n() > 10) %>%
  ungroup()
traits_with_abund_inthigh <- long_abundance_ordered_inthigh %>%
  inner_join(traits_expanded_ordered_inthigh, by = c("species", "row_id"))%>%
  #remove species with less than 4 observations
  group_by(species) %>%
  filter(n() > 10) %>%
  ungroup()

traits_abund <- traits_with_abund %>%
  #remove categorical columns
  select(-Habitat, -Trophic.Level, -Trophic.Niche, -row_id, -longitude, -latitude, -abundance, 
         -Migration, -Primary.Lifestyle)
traits_abund_inthigh <- traits_with_abund_inthigh %>%
  #remove categorical columns
  select(-Habitat, -Trophic.Level, -Trophic.Niche, -row_id, -longitude, -latitude, -abundance, 
         -Migration, -Primary.Lifestyle)

abundance_matrix <- as.matrix(count_predictions[, 3:ncol(count_predictions)]) 
rownames(abundance_matrix) <- paste0("Comm.", seq_len(nrow(abundance_matrix)))

abundance_matrix_inthigh <- as.matrix(count_predictions_inthigh[, 3:ncol(count_predictions_inthigh)]) 
rownames(abundance_matrix_inthigh) <- paste0("Comm.", seq_len(nrow(abundance_matrix_inthigh)))

traits_with_abund <- traits_with_abund %>%
  rename(Hand.Wing.Index = `Hand-Wing.Index`)
traits_with_abund_inthigh <- traits_with_abund_inthigh %>%
  rename(Hand.Wing.Index = `Hand-Wing.Index`)

#----------create dataframes for three different trait combinations (spaces) at regional level------------------
traits_df1 <- data.frame(
  Axis.1 = traits_with_abund$Beak.Length_Culmen,
  Axis.2 = traits_with_abund$Beak.Depth
)
traits_df2 <- data.frame(
  Axis.1 = traits_with_abund$Tail.Length,
  Axis.2 = traits_with_abund$Hand.Wing.Index
)
traits_df3 <- data.frame(
  Axis.1 = traits_with_abund$Habitat.Density,
  Axis.2 = traits_with_abund$Tarsus.Length
)


# Create TPDs and TPDc objects
TPDs <- TPDs(species = traits_with_abund$species, traits_df1)
# Initialize an empty list to store results
TPDc_list <- vector("list", nrow(abundance_matrix))

# Loop through all communities
for (i in 1:nrow(abundance_matrix)) {
  comm_data <- abundance_matrix[i, , drop = FALSE]   # select i-th community
  TPDc_list[[i]] <- TPDc(TPDs = TPDs, sampUnit = comm_data)
}

saveRDS(TPDc_list, file = "TPDc_list1.rds")
TPDc_list <- readRDS("TPDc_list1.rds")


TPDc_list <- vector("list", nrow(abundance_matrix)-2738)


for (i in 2739:nrow(abundance_matrix)) {
  comm_data <- abundance_matrix[i, , drop = FALSE]   # select i-th community
  TPDc_list[[i]] <- TPDc(TPDs = TPDs, sampUnit = comm_data)
}
TPDc_list <- readRDS("TPDc_list2.rds")
TPDc <- TPDc(TPDs = TPDs, sampUnit = abundance_matrix)

# then combine results appropriately

RED <- REND(TPDc = TPDc)
RED_inthigh <- REND(TPDc = TPDc_inthigh)

dissim_d2 <- dissim(TPDs)
dissim_d2_inthigh <- dissim(TPDs_inthigh)

FRed <- redundancy(TPDc = TPDc)
FRed_inthigh <- redundancy(TPDc = TPDc_inthigh)

# Example for a single community
tpd_data <- TPDc$data  # Get the first community's data
grid <- tpd_data$evaluation_grid
prob <- TPDc$TPDc$TPDc$Comm.1

# Build data frame for ggplot

df <- data.frame(grid, prob=prob)
# Mask probability values outside the border so close to zero value are not plotted

df <- df %>%
  mutate(prob = ifelse(prob < 1e-9, NA, prob)) %>%
  mutate(prob = ifelse(
    is.na(prob),
    NA,
    (prob - min(prob, na.rm = TRUE)) /
      (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
  ))

# Normalize grid columns only where prob is not NA

#df<-df%>%drop_na(Axis.1, Axis.2, prob)
# Create a mask for the zero/non-zero areas
df_border <- df %>%
  mutate(has_value = prob > 0.75)  # TRUE if value > 0


# Build FDivergence text block

fr_text <- paste0(
  "Richness = ",
  round(RED$communities$FRichness, 2),
  collapse = "\n"
)

fe_text <- paste0(
  "Evenness = ",
  round(RED$communities$FEvenness, 2),
  collapse = "\n"
)

fd_text <- paste0(
  "Divergence=",
  round(RED$communities$FDivergence, 2),
  collapse = "\n"
)

# Build Redundancy text block
fred_text <- paste(
  paste0(
    "Redundancy=", round(FRed$redundancy, 3)),
  collapse = "\n"
)

# Combine both
annot_text <- paste(fr_text, fe_text, fd_text, fred_text, sep = "\n\n\n\n")