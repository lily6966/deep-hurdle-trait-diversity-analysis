setwd("/Users/liyingnceas/GitHub/Deep-hurdle-population-model-exp/lightning_logs/inthigh")


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
library(exactextractr)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection





#ttt





# Folder path containing raster files
folder_path_historical <- "historical"
folder_path_low <-"prb_change_low"
folder_path_intHigh <-"prb_change_intHigh"

# List all raster files in the folder
#historical_files1 <- list.files(path = folder_path_historical, pattern = ".csv", full.names = TRUE)
low_files <- list.files(path = folder_path_low, pattern = ".csv", full.names = TRUE)
intHigh_files <- list.files(path = folder_path_intHigh, pattern = ".csv", full.names = TRUE)


#--------------------------Presence prob mapping------------------------------------
#create directory to save result
tif_dir <- "occurance_plots_intHigh"
if (!dir.exists(tif_dir)) {
  dir.create(tif_dir)
}
pu <- vect("pu_multi_ccap_small.shp")
# get the buffered checklists for a given year and extract elevation values and calculate median and sd
regions <- pu %>%
  st_as_sf()
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


predictions <- read_csv("filtered_presence_prob_result.csv") #1-12 pevelant species and 13-24 including rare species
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
  ggsave(paste0("occurance_plots_intHigh/", bird_species, "_occurance_probability_historical.pdf"), plot = abundance_change, width = 8, height = 6, units = "in")
}
   


#--------------------------Abundance mapping------------------------------------
#create directory to save result
plot_dir <- "abundance_plots_intHigh"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)

}

predictions <- read_csv("filtered_multiplied_result.csv") #1-12 pevelant species and 13-24 including rare species
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
  

