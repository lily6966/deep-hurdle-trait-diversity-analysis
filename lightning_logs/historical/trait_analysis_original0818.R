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
library(dplyr)
library(tidyr)
library(readxl)
library(exactextractr)
library(TPD)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# === Step 1: Trait matrix ===
# Each row is a species, each column is a trait
count_predictions <-read_csv("multiplied_result.csv") #expectation of abundance
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
  filter(!is.na(abundance) & abundance > 0)

count_filted_inthigh <- count_predictions_inthigh
count_filted_inthigh[ , 3:ncol(count_filted_inthigh)] <- as.data.frame(
  ifelse(
    as.matrix(occurance_binary_inthigh[ , 3:ncol(occurance_binary_inthigh)]) == 1,
    as.matrix(count_predictions_inthigh[ , 3:ncol(count_predictions_inthigh)]),
    NA
  )
)
long_abundance_inthigh <- count_filted_inthigh %>%
  pivot_longer(
    cols = -c(longitude, latitude),
    names_to = "species",
    values_to = "abundance"
  ) %>%
  filter(!is.na(abundance) & abundance > 0)

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
traits_obs_merged <- traits_filtered %>%
  rename(species = Species2) %>%
  inner_join(species_counts, by = "species")

traits_obs_merged_inthigh <- traits_filtered %>%
  rename(species = Species2) %>%
  inner_join(species_counts_inthigh, by = "species")

count_filtered_summary <- count_filtered %>%
  filter(if_any(-c(longitude, latitude), ~ !is.na(.)))%>%
  mutate(across(everything(), ~replace_na(., 0)))%>%
  summarise(across(everything(), ~sum(.))) %>%
  select(-longitude, -latitude)%>%
  select(where(~ any(. != 0)))
count_filtered_summary_inthigh <- count_filted_inthigh %>%
  filter(if_any(-c(longitude, latitude), ~ !is.na(.)))%>%
  #fill NA with 0s
  mutate(across(everything(), ~replace_na(., 0)))%>%
  summarise(across(everything(), ~sum(.))) %>%
  select(-longitude, -latitude)%>%
  select(where(~ any(. != 0)))
traits_obs_filtered <- traits_obs_merged %>%
  filter(species %in% colnames(count_filtered_summary))
traits_obs_filtered_inthigh <- traits_obs_merged_inthigh %>%
  filter(species %in% colnames(count_filtered_summary_inthigh))


# Identify numeric columns (excluding species column)
# Identify trait columns
numeric_cols <- traits_filtered %>%
  select(where(is.numeric)) %>%
  colnames()

categorical_cols <- traits_filtered %>%
  select(-all_of(numeric_cols), -Species2) %>%  # exclude numeric + species
  colnames()
# Expand traits per species
traits_expanded <- do.call(rbind, lapply(1:nrow(traits_obs_filtered), function(i) {
  row <- traits_obs_filtered[i, ]
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

traits_expanded_inthigh <- do.call(rbind, lapply(1:nrow(traits_obs_filtered_inthigh), function(i) {
  row <- traits_obs_filtered_inthigh[i, ]
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

# 1. Arrange both dataframes by species name
long_abundance_ordered <- long_abundance %>%
  arrange(species)

traits_expanded_ordered <- traits_expanded %>%
  arrange(species)

# 2. Create a per-species row index
long_abundance_ordered <- long_abundance_ordered %>%
  group_by(species) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

traits_expanded_ordered <- traits_expanded_ordered %>%
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

long_abundance_ordered_inthigh <- long_abundance_inthigh %>%
  arrange(species)
traits_expanded_ordered_inthigh <- traits_expanded_inthigh %>%
  arrange(species)
long_abundance_ordered_inthigh <- long_abundance_ordered_inthigh %>%
  group_by(species) %>%
  mutate(row_id = row_number()) %>%
  ungroup()
traits_expanded_ordered_inthigh <- traits_expanded_ordered_inthigh %>%
  group_by(species) %>%
  mutate(row_id = row_number()) %>%
  ungroup()
traits_with_abund_inthigh <- long_abundance_ordered_inthigh %>%
  inner_join(traits_expanded_ordered_inthigh, by = c("species", "row_id")) %>%
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
  select(-Habitat, -Trophic.Level, -Trophic.Niche, -row_id, -longitude, -latitude,  -abundance, 
         -Migration, -Primary.Lifestyle)


count_filtered_summary_10 <- count_filtered_summary %>%
  select(intersect(colnames(count_filtered_summary), traits_abund$species)) 
count_filtered_summary_inthigh_10 <- count_filtered_summary_inthigh %>%
  select(intersect(colnames(count_filtered_summary_inthigh), traits_abund_inthigh$species))
# Extract species abundance data
# Combine as two rows by row-binding
abundance_matrix <- as.matrix(count_filtered_summary_10)  # Adjust column range as needed)
# Set column names
#colnames(abundance_matrix) <- names(count_filtered_selected)
abundance_matrix_inthigh <- as.matrix(count_filtered_summary_inthigh_10)
# Assign row names like "Comm.1", "Comm.2", ..., "Comm.5866" required by TPDc
rownames(abundance_matrix) <- paste0("Comm.", seq_len(nrow(abundance_matrix)))
rownames(abundance_matrix_inthigh) <- paste0("Comm.", seq_len(nrow(abundance_matrix_inthigh)))


# 

unique(traits_with_abund$Habitat)
unique(traits_with_abund$Trophic.Level)
unique(traits_with_abund$Trophic.Niche)
unique(traits_with_abund$Primary.Lifestyle)
unique(traits_with_abund$Migration)



# unique(traits_with_abund$Habitat)
# [1] "Wetland"        "Human Modified" "Forest"         "Grassland"      "Riverine"       "Woodland"      
# [7] "Marine"         "Shrubland"      "Coastal"       
# > unique(traits_with_abund$Trophic.Level)
# [1] "Herbivore" "Omnivore"  "Carnivore" "Scavenger"
# > unique(traits_with_abund$Trophic.Niche)
# [1] "Herbivore aquatic"     "Omnivore"              "Granivore"             "Invertivore"          
# [5] "Aquatic predator"      "Scavenger"             "Herbivore terrestrial" "Frugivore"            
# [9] "Vertivore"             "Nectarivore"          
# > unique(traits_with_abund$Primary.Lifestyle)
# [1] "Aquatic"     "Terrestrial" "Insessorial" "Generalist"  "Aerial"     
# > unique(traits_with_abund$Migration)
# [1] "2" "1" "3"

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



traits_df1_inthigh <- data.frame(
  Axis.1 = traits_with_abund_inthigh$Beak.Length_Culmen,
  Axis.2 = traits_with_abund_inthigh$Beak.Depth
)

traits_df2_inthigh <- data.frame(
  Axis.1 = traits_with_abund_inthigh$Tail.Length,
  Axis.2 = traits_with_abund_inthigh$Hand.Wing.Index
)

traits_df3_inthigh <- data.frame(
  Axis.1 = traits_with_abund_inthigh$Habitat.Density,
  Axis.2 = traits_with_abund_inthigh$Tarsus.Length
)



# Create TPDs and TPDc objects
TPDs <- TPDs(species = traits_with_abund$species, traits_df1)
TPDs_inthigh <- TPDs(species = traits_with_abund_inthigh$species, traits_df1_inthigh)

TPDc <- TPDc(TPDs = TPDs, sampUnit = abundance_matrix)
TPDc_inthigh <- TPDc(TPDs = TPDs_inthigh, sampUnit = abundance_matrix_inthigh)


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


# Plot

tpd <- ggplot(df, aes(x = Axis.1, y = Axis.2, fill = prob)) +
  geom_raster() +
  # Add border around non-zero region
  geom_contour(
    data = df_border,
    aes(z = as.numeric(has_value)),
    color = "yellow",
    size = 0.1
  ) +
  #geom_text(data = label_point, aes(x, y, label = "0.25"), inherit.aes = FALSE) +
  scale_fill_gradient(
    low = "lightblue",
    high = "darkblue",
    #limits = c(0.0000001, max(prob)),
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "Beak Length vs Beak Depth",
    x = paste0("Beak Length"),
    y = paste0("Beak Depth")
  ) +
  theme_minimal() +
  #xlim(0, 50) + ylim(0, 20) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),     # axis ticks
    axis.line = element_line(color = "black", linewidth = 0.8) # blac
  ) +
  annotate(
    "text",
    x = min(df$Axis.1, na.rm = TRUE),   # safer to use matching df
    y = max(df$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))

ggsave(paste0("TPC_plot_new/TPDc_beak2.pdf"), plot = tpd, width = 8, height = 6, units = "in", bg = "transparent")


# Example for a single community
tpd_data_inthigh <- TPDc_inthigh$data  # Get the first community's data
grid_inthigh <- tpd_data_inthigh$evaluation_grid
prob_inthigh <- TPDc_inthigh$TPDc$TPDc$Comm.1

# Build data frame for ggplot

df_inthigh <- data.frame(grid_inthigh, prob=prob_inthigh)

df_inthigh <- df_inthigh %>%
  mutate(prob = ifelse(prob < 1e-7, NA, prob)) %>%
  mutate(prob = ifelse(
    is.na(prob),
    NA,
    (prob - min(prob, na.rm = TRUE)) /
      (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
  ))
#df_inthigh <- df_inthigh %>% drop_na(Axis.1, Axis.2, prob)
# Create a mask for the zero/non-zero areas
df_border <- df_inthigh %>%
  mutate(has_value = prob > 0.75)  # TRUE if value > 0

# Build FDivergence text block

fr_text <- paste0(
  "Divergence = ",
  round(RED_inthigh$communities$FRichness, 2),
  collapse = "\n"
)

fe_text <- paste0(
  "Evenness = ",
  round(RED_inthigh$communities$FEvenness, 2),
  collapse = "\n"
)

fd_text <- paste0(
  "Divergence=",
  round(RED_inthigh$communities$FDivergence, 2),
  collapse = "\n"
)

# Build Redundancy text block
fred_text <- paste(
  paste0(
    "Redundancy=", round(FRed_inthigh$redundancy, 3)),
  collapse = "\n"
)

# Combine both
annot_text <- paste(fr_text, fe_text, fd_text, fred_text, sep = "\n\n\n\n")



tpd_inthigh <- ggplot(df_inthigh, aes(x = Axis.1, y = Axis.2, fill = prob)) +
  geom_raster() +
  # Add border around non-zero region
  geom_contour(
    data = df_border,
    aes(z = as.numeric(has_value)),
    color = "red",
    size = 0.4
  ) +
  scale_fill_gradient(
    low = "#FFD6E0",
    high = "#8B004B",
    #limits = c(0.0000001, max(prob)),
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "TPD Plot for Beak Length vs Beak Depth",
    x = paste0("Beak Length (Culmen)"),
    y = paste0("Beak Depth")
  ) +
  theme_minimal() +
  #xlim(0, 50) + ylim(0, 20) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),     # axis ticks
    axis.line = element_line(color = "black", linewidth = 0.8) # blac
  ) +
  annotate(
    "text",
    x = min(df_inthigh$Axis.1, na.rm = TRUE),   # safer to use matching df
    y = max(df_inthigh$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))

#save plot
ggsave(paste0("TPC_plot_new/TPDc_beak1_inthigh.pdf"), plot = tpd_inthigh, width = 8, height = 6, units = "in", bg = "transparent")


# Create TPDs and TPDc objects
TPDs <- TPDs(species = traits_with_abund$species, traits_df2)
TPDs_inthigh <- TPDs(species = traits_with_abund_inthigh$species, traits_df2_inthigh)

TPDc <- TPDc(TPDs = TPDs, sampUnit = abundance_matrix)
TPDc_inthigh <- TPDc(TPDs = TPDs_inthigh, sampUnit = abundance_matrix_inthigh)

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


# Plot

tpd <- ggplot(df, aes(x = grid[, 1], y = grid[, 2], fill = prob)) +
  geom_raster() +
  # Add border around non-zero region
  geom_contour(
    data = df_border,
    aes(z = as.numeric(has_value)),
    color = "yellow",
    size = 0.2
  ) +
  #geom_text(data = label_point, aes(x, y, label = "0.25"), inherit.aes = FALSE) +
  scale_fill_gradient(
    low = "lightblue",
    high = "darkblue",
    #limits = c(0.0000001, max(prob)),
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "Tail Length vs Hand Wing Index",
    x = paste0("Tail Length"),
    y = paste0("Hand Wing Index")
  ) +
  theme_minimal() +
  #xlim(0, 300) + ylim(0, 80) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),     # axis ticks
    axis.line = element_line(color = "black", linewidth = 0.8) # blac
  ) +
  annotate(
    "text",
    x = min(df$Axis.1, na.rm = TRUE),   # safer to use matching df
    y = max(df$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))

ggsave(paste0("TPC_plot_new/TPDc_movement1.pdf"), plot = tpd, width = 8, height = 6, units = "in", bg = "transparent")


# Example for a single community
tpd_data_inthigh <- TPDc_inthigh$data  # Get the first community's data
grid_inthigh <- tpd_data_inthigh$evaluation_grid
prob_inthigh <- TPDc_inthigh$TPDc$TPDc$Comm.1

# Build data frame for ggplot

df_inthigh <- data.frame(grid_inthigh, prob=prob_inthigh)

df_inthigh <- df_inthigh %>%
  mutate(prob = ifelse(prob < 1e-7, NA, prob)) %>%
  mutate(prob = ifelse(
    is.na(prob),
    NA,
    (prob - min(prob, na.rm = TRUE)) /
      (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
  ))

# Create a mask for the zero/non-zero areas
df_border <- df_inthigh %>%
  mutate(has_value = prob > 0.75)  # TRUE if value > 0
# Build FDivergence text block

fr_text <- paste0(
  "Divergence = ",
  round(RED_inthigh$communities$FRichness, 2),
  collapse = "\n"
)

fe_text <- paste0(
  "Evenness = ",
  round(RED_inthigh$communities$FEvenness, 2),
  collapse = "\n"
)

fd_text <- paste0(
  "Divergence=",
  round(RED_inthigh$communities$FDivergence, 2),
  collapse = "\n"
)

# Build Redundancy text block
fred_text <- paste(
  paste0(
    "Redundancy=", round(FRed_inthigh$redundancy, 3)),
  collapse = "\n"
)

# Combine both
annot_text <- paste(fr_text, fe_text, fd_text, fred_text, sep = "\n\n\n\n")



tpd_inthigh <- ggplot(df_inthigh, aes(x = grid_inthigh[, 1], y = grid_inthigh[, 2], fill = prob)) +
  geom_raster() +
  # Add border around non-zero region
  geom_contour(
    data = df_border,
    aes(z = as.numeric(has_value)),
    color = "red",
    size = 0.4
  ) +
  scale_fill_gradient(
    low = "#FFD6E0",
    high = "#8B004B",
    #limits = c(0.0000001, max(prob)),
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "TPD Plot for Tail Length vs Hand Wing Index",
    x = paste0("Tail Length"),
    y = paste0("Hand Wing Index")
  ) +
  theme_minimal() +
  #xlim(0, 300) + ylim(0, 80) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),     # axis ticks
    axis.line = element_line(color = "black", linewidth = 0.8) # blac
  ) +
  annotate(
    "text",
    x = min(df_inthigh$Axis.1, na.rm = TRUE),   # safer to use matching df
    y = max(df_inthigh$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))

#save plot
ggsave(paste0("TPC_plot_new/TPDc_movement1_inthigh.pdf"), plot = tpd_inthigh, width = 8, height = 6, units = "in", bg = "transparent")

# Create TPDs and TPDc objects
TPDs <- TPDs(species = traits_with_abund$species, traits_df3)
TPDs_inthigh <- TPDs(species = traits_with_abund_inthigh$species, traits_df3_inthigh)

TPDc <- TPDc(TPDs = TPDs, sampUnit = abundance_matrix)
TPDc_inthigh <- TPDc(TPDs = TPDs_inthigh, sampUnit = abundance_matrix_inthigh)

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



# Plot

tpd <- ggplot(df, aes(x = grid[, 1], y = grid[, 2], fill = prob)) +
  geom_raster() +
  # Add border around non-zero region
  geom_contour(
    data = df_border,
    aes(z = as.numeric(has_value)),
    color = "yellow",
    size = 0.2
  ) +
  #geom_text(data = label_point, aes(x, y, label = "0.25"), inherit.aes = FALSE) +
  scale_fill_gradient(
    low = "lightblue",
    high = "darkblue",
    #limits = c(0.0000001, max(prob)),
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "Habitat Density vs Tarsus Length",
    x = paste0("Habitat Density"),
    y = paste0("Tarsus Length")
  ) +
  theme_minimal() +
  #xlim(0, 3) + ylim(0, 100) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),     # axis ticks
    axis.line = element_line(color = "black", linewidth = 0.8) # blac
  ) +
  annotate(
    "text",
    x = min(df$Axis.1, na.rm = TRUE),   # safer to use matching df
    y = max(df$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))

ggsave(paste0("TPC_plot_new/TPDc_size1.pdf"), plot = tpd, width = 8, height = 6, units = "in", bg = "transparent")


# Example for a single community
tpd_data_inthigh <- TPDc_inthigh$data  # Get the first community's data
grid_inthigh <- tpd_data_inthigh$evaluation_grid
prob_inthigh <- TPDc_inthigh$TPDc$TPDc$Comm.1

# Build data frame for ggplot

df_inthigh <- data.frame(grid_inthigh, prob=prob_inthigh)

df_inthigh <- df_inthigh %>%
  mutate(prob = ifelse(prob < 1e-7, NA, prob)) %>%
  mutate(prob = ifelse(
    is.na(prob),
    NA,
    (prob - min(prob, na.rm = TRUE)) /
      (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
  ))

# Create a mask for the zero/non-zero areas
df_border <- df_inthigh %>%
  mutate(has_value = prob > 0.75)  # TRUE if value > 0
# Build FDivergence text block

fr_text <- paste0(
  "Richness = ",
  round(RED_inthigh$communities$FRichness, 2),
  collapse = "\n"
)

fe_text <- paste0(
  "Evenness = ",
  round(RED_inthigh$communities$FEvenness, 2),
  collapse = "\n"
)

fd_text <- paste0(
  "Divergence=",
  round(RED_inthigh$communities$FDivergence, 2),
  collapse = "\n"
)

# Build Redundancy text block
fred_text <- paste(
  paste0(
    "Redundancy=", round(FRed_inthigh$redundancy, 3)),
  collapse = "\n"
)

# Combine both
annot_text <- paste(fr_text, fe_text, fd_text, fred_text, sep = "\n\n\n\n")



tpd_inthigh <- ggplot(df_inthigh, aes(x = grid_inthigh[, 1], y = grid_inthigh[, 2], fill = prob)) +
  geom_raster() +
  # Add border around non-zero region
  geom_contour(
    data = df_border,
    aes(z = as.numeric(has_value)),
    color = "red",
    size = 0.4
  ) +
  scale_fill_gradient(
    low = "#FFD6E0",
    high = "#8B004B",
    #limits = c(0.0000001, max(prob)),
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "TPD Plot for Size vs Habitat Density",
    x = paste0("Habitat Density"),
    y = paste0("Tarsus Length")
  ) +
  theme_minimal() +
  #xlim(0, 3) + ylim(0, 100) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),     # axis ticks
    axis.line = element_line(color = "black", linewidth = 0.8) # blac
  ) +
  annotate(
    "text",
    x = min(df_inthigh$Axis.1, na.rm = TRUE),   # safer to use matching df
    y = max(df_inthigh$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))

#save plot
ggsave(paste0("TPC_plot_new/TPDc_size1_inthigh.pdf"), plot = tpd_inthigh, width = 8, height = 6, units = "in", bg = "transparent")

#by trait categories
# Step 3: Plot

for (cat_col in categorical_cols) {

  # Get the unique values of the current column (excluding NA)
  unique_vals <- na.omit(unique(traits_with_abund[[cat_col]]))

  for (val in unique_vals) {
    # Filter the data for the current category value
    filtered_data <- traits_with_abund %>%
      filter(!!sym(cat_col) == val)


      traits_abund_element <- traits_with_abund %>%
        filter(!!sym(cat_col) == val) %>%
        #remove categorical columns
        select(-longitude, -latitude,-Habitat, -Trophic.Level, -Trophic.Niche,  -abundance,
               -Migration, -Primary.Lifestyle, -row_id)


      traits_abund_element_inthigh <- traits_with_abund_inthigh %>%
        filter(!!sym(cat_col) == val) %>%
        #remove categorical columns
        select(-longitude, -latitude,-Habitat, -Trophic.Level, -Trophic.Niche,  -abundance,
                 -Migration, -Primary.Lifestyle, -row_id)

      # Skip if too big
      count_filtered_selected <- count_filtered %>%
        select(all_of(traits_abund_element$species)) %>%
        #fill NA with 0s
        mutate(across(everything(), ~replace_na(., 0)))%>%
        summarise(across(everything(), ~sum(.)))
      count_filtered_selected_inthigh <- count_filted_inthigh %>%
        filter(if_any(-c(longitude, latitude), ~ !is.na(.)))%>%
        select(all_of(traits_abund_element_inthigh$species)) %>%
        #fill NA with 0s
        mutate(across(everything(), ~replace_na(., 0)))%>%
        summarise(across(everything(), ~sum(.)))


      # Extract species abundance data
      # Combine as two rows by row-binding
      abundance_matrix <- as.matrix(count_filtered_selected)  # Adjust column range as needed)
      # Set column names
      #colnames(abundance_matrix) <- names(count_filtered_selected)
      abundance_matrix_inthigh <- as.matrix(count_filtered_selected_inthigh)
      # Assign row names like "Comm.1", "Comm.2", ..., "Comm.5866" required by TPDc
      rownames(abundance_matrix) <- paste0("Comm.", seq_len(nrow(abundance_matrix)))
      rownames(abundance_matrix_inthigh) <- paste0("Comm.", seq_len(nrow(abundance_matrix_inthigh)))




        traits_df1 <- data.frame(
          Axis.1 = traits_abund_element$Beak.Length_Culmen,
          Axis.2 = traits_abund_element$Beak.Depth
        )
        traits_df2 <- data.frame(
          Axis.1 = traits_abund_element$Tail.Length,
          Axis.2 = traits_abund_element$Hand.Wing.Index
        )
        traits_df3 <- data.frame(
          Axis.1 = traits_abund_element$Habitat.Density,
          Axis.2 = traits_abund_element$Tarsus.Length
        )

        traits_df1_inthigh <- data.frame(
          Axis.1 = traits_abund_element_inthigh$Beak.Length_Culmen,
          Axis.2 = traits_abund_element_inthigh$Beak.Depth
        )
        traits_df2_inthigh <- data.frame(
          Axis.1 = traits_abund_element_inthigh$Tail.Length,
          Axis.2 = traits_abund_element_inthigh$Hand.Wing.Index
        )
        traits_df3_inthigh <- data.frame(
          Axis.1 = traits_abund_element_inthigh$Habitat.Density,
          Axis.2 = traits_abund_element_inthigh$Tarsus.Length
        )


    # Create TPDs and TPDc objects
    TPDs_element <- TPDs(species = traits_abund_element$species, traits_df1)
    TPDs_element_inthigh <- TPDs(species = traits_abund_element_inthigh$species, traits_df1_inthigh)

    TPDc_element <- TPDc(TPDs = TPDs_element, sampUnit = abundance_matrix)
    TPDc_element_inthigh <- TPDc(TPDs = TPDs_element_inthigh, sampUnit = abundance_matrix_inthigh)

    # Example for a single community
    tpd_data <- TPDc_element$data  # Get the first community's data
    grid <- tpd_data$evaluation_grid
    prob <- TPDc_element$TPDc$TPDc$Comm.1

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

    # Create a mask for the zero/non-zero areas
    df_border <- df %>%
      mutate(has_value = prob > 0.75)  # TRUE if value > 0

    # Plot

    tpd_element <- ggplot(df, aes(x = grid[, 1], y = grid[, 2], fill = prob)) +
      geom_raster() +
      # Add border around non-zero region
      geom_contour(
        data = df_border,
        aes(z = as.numeric(has_value)),
        color = "yellow",
        size = 0.2
      ) +
      #geom_text(data = label_point, aes(x, y, label = "0.25"), inherit.aes = FALSE) +
      scale_fill_gradient(
        low = "lightblue",
        high = "darkblue",
        #limits = c(0.0000001, max(prob)),
        name = "TPD Value",
        na.value = "transparent"
      ) +
      labs(
        title = "Beak Length vs Beak Depth TPD",
        x = paste0("Beak Length"),
        y = paste0("Beak Depth")
      ) +
      theme_minimal() +
      #xlim(0, 50) + ylim(0, 20) +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),  # remove major gridlines
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),     # axis ticks
        axis.line = element_line(color = "black", linewidth = 0.8) # blac
      ) +
      guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))

    ggsave(paste0("TPC_plot_new/TPDc_beak_", val,".pdf"), plot = tpd_element, width = 8, height = 6, units = "in", bg = "transparent")


    # Example for a single community
    tpd_data_inthigh <- TPDc_element_inthigh$data  # Get the first community's data
    grid_inthigh <- tpd_data_inthigh$evaluation_grid
    prob_inthigh <- TPDc_element_inthigh$TPDc$TPDc$Comm.1

    # Build data frame for ggplot

    df_inthigh <- data.frame(grid_inthigh, prob=prob_inthigh)

    df_inthigh <- df_inthigh %>%
      mutate(prob = ifelse(prob < 1e-7, NA, prob)) %>%
      mutate(prob = ifelse(
        is.na(prob),
        NA,
        (prob - min(prob, na.rm = TRUE)) /
          (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
      ))

    # Create a mask for the zero/non-zero areas
    df_border <- df_inthigh %>%
      mutate(has_value = prob > 0.75)  # TRUE if value > 0


    tpd_element_inthigh <- ggplot(df_inthigh, aes(x = grid_inthigh[, 1], y = grid_inthigh[, 2], fill = prob)) +
      geom_raster() +
      # Add border around non-zero region
      geom_contour(
        data = df_border,
        aes(z = as.numeric(has_value)),
        color = "red",
        size = 0.4
      ) +
      scale_fill_gradient(
        low = "#FFD6E0",
        high = "#8B004B",
        #limits = c(0.0000001, max(prob)),
        name = "TPD Value",
        na.value = "transparent"
      ) +
      labs(
        title = "TPD Plot for Beak Length vs Beak Depth",
        x = paste0("Beak Length"),
        y = paste0("Beak Depth")
      ) +
      theme_minimal() +
      #xlim(0, 50) + ylim(0, 20) +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),  # remove major gridlines
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),     # axis ticks
        axis.line = element_line(color = "black", linewidth = 0.8) # blac
      ) +
      guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))

    #save plot
    ggsave(paste0("TPC_plot_new/TPDc_beak_", val, "_inthigh.pdf"), plot = tpd_element_inthigh, width = 8, height = 6, units = "in", bg = "transparent")


    # Create TPDs and TPDc objects
    TPDs_element <- TPDs(species = traits_abund_element$species, traits_df3)
    TPDs_element_inthigh <- TPDs(species = traits_abund_element_inthigh$species, traits_df3_inthigh)
    
    TPDc_element <- TPDc(TPDs = TPDs_element, sampUnit = abundance_matrix)
    TPDc_element_inthigh <- TPDc(TPDs = TPDs_element_inthigh, sampUnit = abundance_matrix_inthigh)
    
    # Example for a single community
    tpd_data <- TPDc_element$data  # Get the first community's data
    grid <- tpd_data$evaluation_grid
    prob <- TPDc_element$TPDc$TPDc$Comm.1
    
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
    
    # Create a mask for the zero/non-zero areas
    df_border <- df %>%
      mutate(has_value = prob > 0.75)  # TRUE if value > 0
    
    # Plot
    
    tpd_element <- ggplot(df, aes(x = grid[, 1], y = grid[, 2], fill = prob)) +
      geom_raster() +
      # Add border around non-zero region
      geom_contour(
        data = df_border,
        aes(z = as.numeric(has_value)),
        color = "yellow",
        size = 0.2
      ) +
      #geom_text(data = label_point, aes(x, y, label = "0.25"), inherit.aes = FALSE) +
      scale_fill_gradient(
        low = "lightblue",
        high = "darkblue",
        #limits = c(0.0000001, max(prob)),
        name = "TPD Value",
        na.value = "transparent"
      ) +
      labs(
        title = "Habitat Density vs Tarsus Length TPD",
        x = paste0("Habitat Density"),
        y = paste0("Tarsus Length")
      ) +
      theme_minimal() +
      #xlim(0, 3) + ylim(0, 100) +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),  # remove major gridlines
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(color = "black"),     # axis ticks
        axis.line = element_line(color = "black", linewidth = 0.8) # blac
      ) +
      guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))
    
    ggsave(paste0("TPC_plot_new/TPDc_size_", val,".pdf"), plot = tpd_element, width = 8, height = 6, units = "in", bg = "transparent")
    
    
    # Example for a single community
    tpd_data_inthigh <- TPDc_element_inthigh$data  # Get the first community's data
    grid_inthigh <- tpd_data_inthigh$evaluation_grid
    prob_inthigh <- TPDc_element_inthigh$TPDc$TPDc$Comm.1
    
    # Build data frame for ggplot
    
    df_inthigh <- data.frame(grid_inthigh, prob=prob_inthigh)
    
    df_inthigh <- df_inthigh %>%
      mutate(prob = ifelse(prob < 1e-7, NA, prob)) %>%
      mutate(prob = ifelse(
        is.na(prob), 
        NA, 
        (prob - min(prob, na.rm = TRUE)) /
          (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
      ))
    
    # Create a mask for the zero/non-zero areas
    df_border <- df_inthigh %>%
      mutate(has_value = prob > 0.75)  # TRUE if value > 0
    
    
    tpd_element_inthigh <- ggplot(df_inthigh, aes(x = grid_inthigh[, 1], y = grid_inthigh[, 2], fill = prob)) +
      geom_raster() +
      # Add border around non-zero region
      geom_contour(
        data = df_border,
        aes(z = as.numeric(has_value)),
        color = "red",
        size = 0.4
      ) +
      scale_fill_gradient(
        low = "#FFD6E0",
        high = "#8B004B",
        #limits = c(0.0000001, max(prob)),
        name = "TPD Value",
        na.value = "transparent"
      ) +
      labs(
        title = "TPD Plot for Mass vs Habitat Density",
        x = paste0("Mass"),
        y = paste0("Habitat Density")
      ) +
      theme_minimal() +
      #xlim(0, 3) + ylim(0, 100) +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),  # remove major gridlines
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(color = "black"),     # axis ticks
        axis.line = element_line(color = "black", linewidth = 0.8) # blac
      ) +
      guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))
    
    #save plot
    ggsave(paste0("TPC_plot_new/TPDc_size_", val, "_inthigh.pdf"), plot = tpd_element_inthigh, width = 8, height = 6, units = "in", bg = "transparent")
    
    
  } 
}

for (cat_col in categorical_cols) {
  
  # Get the unique values of the current column (excluding NA)
  unique_vals <- na.omit(unique(traits_with_abund[[cat_col]]))
  
  for (val in unique_vals) {
    # Filter the data for the current category value
    filtered_data <- traits_with_abund %>%
      filter(!!sym(cat_col) == val)
    traits_abund_element <- traits_with_abund %>%
      filter(!!sym(cat_col) == val) %>%
      #remove categorical columns
      select(-Habitat, -Trophic.Level, -Trophic.Niche, -row_id)
    traits_abund_element_inthigh <- traits_with_abund_inthigh %>%
      filter(!!sym(cat_col) == val) %>%
      #remove categorical columns
      select(-Habitat, -Trophic.Level, -Trophic.Niche, -row_id)
    # Split by species
    traits_split <- split(traits_abund_element, traits_abund_element$species)
    
    
    gower_dist <- lapply(traits_split, function(df) {
      df_traits <- df %>% 
        select(-longitude, -latitude, -species, -abundance, 
               -Migration, -Primary.Lifestyle) %>%
        mutate(across(where(is.character), as.factor)) %>%
        mutate(across(where(is.logical), as.numeric)) %>%
        mutate(across(everything(), ~ {
          # Check if variable has exactly 2 unique values (excluding NAs)
          unique_vals <- unique(na.omit(.))
          if (length(unique_vals) == 2 && !is.numeric(.)) {
            as.numeric(as.character(.))
          } else {
            .
          }
        }))
      
      if (nrow(df_traits) > 1) {
        tryCatch({
          gowdis(df_traits)
        }, error = function(e) {
          message("Error for species: ", unique(df$species), " - ", e$message)
          return(NULL)
        })
      } else {
        NULL
      }
    })
    
    library(ape)
    
    # gower_dist is a named list of distance matrices, one per species
    
    pcoa_list <- lapply(names(gower_dist), function(species_name) {
      dist_mat <- gower_dist[[species_name]]
      
      if (is.null(dist_mat)) {
        return(NULL)  # skip species with no distance matrix
      }
      
      # Run PCoA
      pcoa_res <- pcoa(dist_mat)
      
      # Extract first two PCoA axes
      scores <- as.data.frame(pcoa_res$vectors[, 1:2])
      colnames(scores) <- c("PCoA1", "PCoA2")
      
      # Add species info and individual IDs
      scores$species <- species_name
      scores$individual_id <- rownames(dist_mat)
      
      return(scores)
    })
    
    # Combine all into one data frame, removing NULLs
    all_pcoa_scores <- bind_rows(pcoa_list)
    
    pcoa_df <- as.data.frame(all_pcoa_scores)
    
    # Axis	Relative Eigenvalue	Cumulative Variance
    # 1	0.4407	0.4407
    # 2	0.2452	0.6860
    # 3	0.1139	0.7999
    # 4	0.0752	0.8751
    pcoa_2traits_df <- pcoa_df %>% select(PCoA1, PCoA2)
    colnames(pcoa_2traits_df) <- paste0("Axis.", 1:2)
    # Split by species
    traits_split_inthigh <- split(traits_abund_element_inthigh, traits_abund_element_inthigh$species)
    
    
    gower_dist_inthigh <- lapply(traits_split_inthigh, function(df) {
      df_traits <- df %>% 
        select(-longitude, -latitude, -species, -abundance, 
               -Migration, -Primary.Lifestyle) %>%
        mutate(across(where(is.character), as.factor)) %>%
        mutate(across(where(is.logical), as.numeric)) %>%
        mutate(across(everything(), ~ {
          # Check if variable has exactly 2 unique values (excluding NAs)
          unique_vals <- unique(na.omit(.))
          if (length(unique_vals) == 2 && !is.numeric(.)) {
            as.numeric(as.character(.))
          } else {
            .
          }
        }))
      
      if (nrow(df_traits) > 1) {
        tryCatch({
          gowdis(df_traits)
        }, error = function(e) {
          message("Error for species: ", unique(df$species), " - ", e$message)
          return(NULL)
        })
      } else {
        NULL
      }
    })
    
    library(ape)
    
    # gower_dist is a named list of distance matrices, one per species
    
    pcoa_list_inthigh <- lapply(names(gower_dist_inthigh), function(species_name) {
      dist_mat <- gower_dist_inthigh[[species_name]]
      
      if (is.null(dist_mat)) {
        return(NULL)  # skip species with no distance matrix
      }
      
      # Run PCoA
      pcoa_res <- pcoa(dist_mat)
      
      # Extract first two PCoA axes
      scores <- as.data.frame(pcoa_res$vectors[, 1:2])
      colnames(scores) <- c("PCoA1", "PCoA2")
      
      # Add species info and individual IDs
      scores$species <- species_name
      scores$individual_id <- rownames(dist_mat)
      
      return(scores)
    })
    
    # Combine all into one data frame, removing NULLs
    all_pcoa_scores_inthigh <- bind_rows(pcoa_list_inthigh)
    
    pcoa_df_inthigh <- as.data.frame(all_pcoa_scores_inthigh)
    
    # Axis	Relative Eigenvalue	Cumulative Variance
    # 1	0.4407	0.4407
    # 2	0.2452	0.6860
    # 3	0.1139	0.7999
    # 4	0.0752	0.8751
    pcoa_2traits_df_inthigh <- pcoa_df_inthigh %>% select(PCoA1, PCoA2)
    
    colnames(pcoa_2traits_df_inthigh) <- paste0("Axis.", 1:2)
    TPDs_element <- TPDs(species = traits_abund_element$species, pcoa_2traits_df)
    TPDs_element_inthigh <- TPDs(species = traits_abund_element_inthigh$species, pcoa_2traits_df_inthigh)
    
    
    count_filtered_selected <- count_filtered %>%
      select(all_of(traits_abund_element$species)) %>%
      #fill NA with 0s
      mutate(across(everything(), ~replace_na(., 0)))%>%
      summarise(across(everything(), ~sum(.)))
    count_filtered_selected_inthigh <- count_filted_inthigh %>%
      filter(if_any(-c(longitude, latitude), ~ !is.na(.)))%>%
      select(all_of(traits_abund_element_inthigh$species)) %>%
      #fill NA with 0s
      mutate(across(everything(), ~replace_na(., 0)))%>%
      summarise(across(everything(), ~sum(.)))
    
    # Extract species abundance data
    # Combine as two rows by row-binding
    abundance_matrix <- as.matrix(count_filtered_selected)  # Adjust column range as needed)
    # Set column names
    #colnames(abundance_matrix) <- names(count_filtered_selected)
    abundance_matrix_inthigh <- as.matrix(count_filtered_selected_inthigh)
    # Assign row names like "Comm.1", "Comm.2", ..., "Comm.5866" required by TPDc
    rownames(abundance_matrix) <- paste0("Comm.", seq_len(nrow(abundance_matrix)))
    rownames(abundance_matrix_inthigh) <- paste0("Comm.", seq_len(nrow(abundance_matrix_inthigh)))
    TPDc_element <- TPDc(TPDs = TPDs_element, sampUnit = abundance_matrix)
    TPDc_element_inthigh <- TPDc(TPDs = TPDs_element_inthigh, sampUnit = abundance_matrix_inthigh)
    
    # Example for a single community
    tpd_data <- TPDc_element$data  # Get the first community's data
    grid <- tpd_data$evaluation_grid
    prob <- TPDc_element$TPDc$TPDc$Comm.1
    
    # Build data frame for ggplot
    
    df <- data.frame(grid, prob=prob)
    # Mask probability values outside the border so close to zero value are not plotted
    df <- df %>%
      mutate(prob = ifelse(prob < 0.0000001, NA, prob))
    
    # Create a mask for the zero/non-zero areas
    df_border <- df %>%
      mutate(has_value = prob > 0.5)  # TRUE if value > 0
    
    tpd_element <- ggplot(df, aes(x = grid[, 1], y = grid[, 2], fill = prob)) +
      geom_raster() +
      # Add border around non-zero region
      geom_contour(
        data = df_border,
        aes(z = as.numeric(has_value)),
        color = "black",
        size = 0.4
      ) +
      scale_fill_gradient(
        low = "lightblue",
        high = "darkblue",
        limits = c(0.0000001, max(prob)),
        name = "TPD Value",
        na.value = "transparent"
      ) +
      labs(
        title = "TPD Plot",
        x = "Axes 1",
        y = "Axes 2"
      ) +
      theme_minimal() +
      coord_equal() +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),  # remove major gridlines
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(color = "black"),     # axis ticks
        axis.line = element_line(color = "black", linewidth = 0.8) # blac
      ) +
      guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))
    
    ggsave(paste0("TPC_plot_new/TPDc_", val,".pdf"), plot = tpd_element, width = 8, height = 6, units = "in", bg = "transparent")
    
    
    # Example for a single community
    tpd_data_inthigh <- TPDc_element_inthigh$data  # Get the first community's data
    grid_inthigh <- tpd_data_inthigh$evaluation_grid
    prob_inthigh <- TPDc_element_inthigh$TPDc$TPDc$Comm.1
    
    # Build data frame for ggplot
    
    df_inthigh <- data.frame(grid_inthigh, prob=prob_inthigh)
    
    
    # Mask probability values outside the border so close to zero value are not plotted
    df_inthigh <- df_inthigh %>%
      mutate(prob = ifelse(prob < 0.0000001, NA, prob))
    
    tpd_element_inthigh <- ggplot(df_inthigh, aes(x = grid_inthigh[, 1], y = grid_inthigh[, 2], fill = prob)) +
      geom_raster() +
      # Add border around non-zero region
      # geom_contour(
      #   data = df_border,
      #   aes(z = as.numeric(has_value)),
      #   color = "black",
      #   size = 0.4
      # ) +
      scale_fill_gradient(
        low = "#FFD6E0",
        high = "#8B004B",
        limits = c(0.0000001, max(prob)),
        name = "TPD Value",
        na.value = "transparent"
      ) +
      labs(
        title = "TPD Plot",
        x = "Axes 1",
        y = "Axes 2"
      ) +
      theme_minimal() +
      coord_equal() +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),  # remove major gridlines
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(color = "black"),     # axis ticks
        axis.line = element_line(color = "black", linewidth = 0.8) # blac
      ) +
      guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))
    
    #save plot
    ggsave(paste0("TPC_plot_new/TPDc_", val, "_inthigh.pdf"), plot = tpd_element_inthigh, width = 8, height = 6, units = "in", bg = "transparent")
    
    
    
  } 
}

library(terra)
library(sf)

# Read shapefile (.shp)
shp_apa <- st_read("RegionBoundaries/APA_boundary.shp")%>%
  # project to the cdl projection
  st_transform(crs = 4326)
shp_gb <- st_read("RegionBoundaries/GB_boundary.shp")%>%
  # project to the cdl projection
  st_transform(crs = 4326)
shp_wb <- st_read("RegionBoundaries/WB_boundary.shp")%>%
  # project to the cdl projection
  st_transform(crs = 4326)

#remove rows only with NAs
count_filtered_clean <- count_filtered %>%
  filter(if_any(-c(longitude, latitude), ~ !is.na(.)))

# Keep only longitude, latitude, and species columns in `species_keep`
count_filtered_clean <- count_filtered_clean %>%
  
  #fill NA with 0s
  mutate(across(everything(), ~replace_na(., 0)))

count_filtered_clean_inthigh <- count_filted_inthigh %>%
  filter(if_any(-c(longitude, latitude), ~ !is.na(.)))%>%
  #fill NA with 0s
  mutate(across(everything(), ~replace_na(., 0)))

# 1. Convert your tibble to an sf POINT object using lon/lat
points_sf <- count_filtered_clean %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

points_sf_inthigh <- count_filtered_clean_inthigh %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# 2. Filter points that fall inside the polygon
# points_in_apa <- points_sf[st_within(points_sf, shp_apa, sparse = FALSE), ]
# points_in_gb <- points_sf[st_within(points_sf, shp_gb, sparse = FALSE), ]
# points_in_wb <- points_sf[st_within(points_sf, shp_wb, sparse = FALSE), ]
points_in_any <- points_sf[
  st_within(points_sf, shp_apa, sparse = FALSE) |
    st_within(points_sf, shp_gb, sparse = FALSE) |
    st_within(points_sf, shp_wb, sparse = FALSE),
]


#points_in_apa_inthigh <- points_sf_inthigh[st_within(points_sf_inthigh, shp_apa, sparse = FALSE), ]
#points_in_gb_inthigh <- points_sf_inthigh[st_within(points_sf_inthigh, shp_gb, sparse = FALSE), ]
#points_in_wb_inthigh <- points_sf_inthigh[st_within(points_sf_inthigh, shp_wb, sparse = FALSE), ]
points_in_any_inthigh <- points_sf_inthigh[
  st_within(points_sf_inthigh, shp_apa, sparse = FALSE) |
    st_within(points_sf_inthigh, shp_gb, sparse = FALSE) |
    st_within(points_sf_inthigh, shp_wb, sparse = FALSE),
]

sf_in_any <- points_in_any %>%
  mutate(longitude = st_coordinates(points_in_any)[, "X"],
         latitude = st_coordinates(points_in_any)[, "Y"]) %>%
  st_drop_geometry()

sf_in_any_inthigh <- points_in_any_inthigh %>%
  mutate(longitude = st_coordinates(points_in_any_inthigh)[, "X"],
         latitude = st_coordinates(points_in_any_inthigh)[, "Y"]) %>%
  st_drop_geometry()





# Pivot to long format: one row per site-species pair filtered by species occurance
long_abundance_comms <- sf_in_any %>%
  
  pivot_longer(
    cols = -c(longitude, latitude),
    names_to = "species",
    values_to = "abundance"
  ) %>% filter (abundance > 0)
long_abundance_comms_inthigh <- sf_in_any_inthigh %>%
  pivot_longer(
    cols = -c(longitude, latitude),
    names_to = "species",
    values_to = "abundance"
  )%>% filter (abundance > 0)
species_counts_comms <- long_abundance_comms %>%
  count(species, name = "n_obs")
species_counts_comms_inthigh <- long_abundance_comms_inthigh %>%
  count(species, name = "n_obs")

#merge trait means with observation counts
trait_obs_merged_comms <- traits_filtered %>%
  rename(species = Species2) %>%
  inner_join(species_counts_comms, by = "species")

trait_obs_merged_comms_inthigh <- traits_filtered %>%
  rename(species = Species2) %>%
  inner_join(species_counts_inthigh, by = "species")

numeric_cols <- traits_filtered %>%
  select(where(is.numeric)) %>%
  colnames()
categorical_cols <- traits_filtered %>%
  select(where(is_character)) %>%
  colnames()


# auto-detect numeric columns and exclude the 'n_obs' counter
numeric_cols <- names(trait_obs_merged_comms_inthigh)[vapply(trait_obs_merged_comms_inthigh, is.numeric, TRUE)]
numeric_cols <- setdiff(numeric_cols, "n_obs")

# auto-detect categorical columns (non-numeric, excluding species)
categorical_cols <- setdiff(names(trait_obs_merged_comms_inthigh)[!vapply(trait_obs_merged_comms_inthigh, is.numeric, TRUE)],
                            c("species"))

# Expand traits per species
traits_expanded_2_comms <- do.call(rbind, lapply(1:nrow(trait_obs_merged_comms), function(i) {
  row <- trait_obs_merged_comms[i, ]
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

traits_expanded_2_comms_inthigh <- do.call(rbind, lapply(1:nrow(trait_obs_merged_comms_inthigh), function(i) {
  row <- trait_obs_merged_comms_inthigh[i, ]
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




# 1. Arrange both dataframes by species name
long_abundance_ordered_comms <- long_abundance_comms %>%
  arrange(species)

traits_expanded_ordered_comms <- traits_expanded_2_comms %>%
  arrange(species)


long_abundance_ordered_comms_inthigh <- long_abundance_inthigh %>%
  arrange(species)

traits_expanded_ordered_comms_inthigh <- traits_expanded_2_comms_inthigh %>%
  arrange(species)

# 2. Create a per-species row index
long_abundance_ordered_comms <- long_abundance_ordered_comms %>%
  group_by(species) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

traits_expanded_ordered_comms <- traits_expanded_ordered_comms %>%
  group_by(species) %>%
  mutate(row_id = row_number()) %>%
  ungroup()


long_abundance_ordered_comms_inthigh <- long_abundance_ordered_comms_inthigh %>%
  group_by(species) %>%
  mutate(row_id = row_number()) %>%
  ungroup()
traits_expanded_ordered_comms_inthigh <- traits_expanded_ordered_comms_inthigh %>%
  group_by(species) %>%
  mutate(row_id = row_number()) %>%
  ungroup()
# 3. Join by species and row_id (ensures correct alignment)
traits_with_abund_2comms <- long_abundance_ordered_comms %>%
  inner_join(traits_expanded_ordered_comms, by = c("species", "row_id"))%>%
  #remove species with less than 4 observations
  group_by(species) %>%
  filter(n() > 10) %>%
  ungroup()

traits_with_abund_2comms_inthigh <- long_abundance_ordered_comms_inthigh %>%
  inner_join(traits_expanded_ordered_comms_inthigh, by = c("species", "row_id"))%>%
  #remove species with less than 4 observations
  group_by(species) %>%
  filter(n() > 10) %>%
  ungroup()

traits_with_abund_2comms <- traits_with_abund_2comms %>%
  rename(Hand.Wing.Index = 'Hand-Wing.Index') 
traits_with_abund_2comms_inthigh <- traits_with_abund_2comms_inthigh %>%
  rename(Hand.Wing.Index = 'Hand-Wing.Index')

#create foraging trait dataframes
traits_df1_comms <- data.frame(
  Axis.1 =traits_with_abund_2comms$Beak.Length_Culmen,
  Axis.2 = traits_with_abund_2comms$Beak.Depth
)

traits_df1_comms_inthigh <- data.frame(
  Axis.1 =traits_with_abund_2comms_inthigh$Beak.Length_Culmen,
  Axis.2 = traits_with_abund_2comms_inthigh$Beak.Depth
)

#create movement trait dataframes
traits_df2_comms <- data.frame(
  Axis.1 = traits_with_abund_2comms$Tail.Length,
  Axis.2 = traits_with_abund_2comms$Hand.Wing.Index
)

traits_df2_comms_inthigh <- data.frame(
  Axis.1 = traits_with_abund_2comms_inthigh$Tail.Length,
  Axis.2 = traits_with_abund_2comms_inthigh$Hand.Wing.Index
)

#create social behavior trait dataframes
traits_df3_comms <- data.frame(
  Axis.1 = traits_with_abund_2comms$Habitat.Density,
  Axis.2 = traits_with_abund_2comms$Tarsus.Length
)

traits_df3_comms_inthigh <- data.frame(
  Axis.1 = traits_with_abund_2comms_inthigh$Habitat.Density,
  Axis.2 = traits_with_abund_2comms_inthigh$Tarsus.Length
)


# beak plots
TPDs_d2 <- TPDs(species = traits_with_abund_2comms$species, traits_df1_comms)
RED_d2 <- REND(TPDs = TPDs_d2)



TPDs_d2_inthigh <- TPDs(species = traits_with_abund_2comms_inthigh$species, traits_df1_comms_inthigh)
RED_d2_inthigh <- REND(TPDs = TPDs_d2_inthigh)

# -----Create community abundance matrix for two communities for all trait spaces------------------------
df_two_comms <- bind_rows(
  sf_in_any %>% mutate(region = "Historical"),
  sf_in_any_inthigh %>% mutate(region = "Sea Level Rise")
) %>%
  mutate(region = factor(region, levels = c("Historical", "Sea Level Rise"))) %>%
  select(-longitude, -latitude) %>%
  group_by(region) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  select(any_of(unique(traits_with_abund_2comms_inthigh$species)))
df_two_matrix<- as.matrix(df_two_comms)
rownames(df_two_matrix) <- paste0("Comm.", seq_len(nrow(df_two_matrix)))

sampUnit <- df_two_matrix[1, , drop = FALSE] # 1-row dataframe, species columns only
sampUnit_inthigh <- df_two_matrix[2, , drop = FALSE]
#-------------------------------------------------------------------------------------------------------
#-------Start TPDc and REND calculations for beak trait space reserve level------------------------------
TPDc_divergence_d2 <- TPDc(TPDs = TPDs_d2, sampUnit = sampUnit)
RED_comm_d2 <- REND(TPDc = TPDc_divergence_d2)

TPDc_divergence_d2_inthigh <- TPDc(TPDs = TPDs_d2_inthigh, sampUnit = sampUnit_inthigh)
RED_comm_d2_inthigh <- REND(TPDc = TPDc_divergence_d2_inthigh)

dissim_d2 <- dissim(TPDs_d2)
dissim_d2_inthigh <- dissim(TPDs_d2_inthigh)

FRed_d2 <- redundancy(TPDc = TPDc_divergence_d2)
FRed_d2_inthigh <- redundancy(TPDc = TPDc_divergence_d2_inthigh)


# Example for a single community
tpd_data_ <- TPDc_divergence_d2$data  # Get the data
grid <- tpd_data_$evaluation_grid
prob <- TPDc_divergence_d2$TPDc$TPDc$Comm.1

# Build data frame for ggplot
df_comms <- data.frame(grid, prob = prob)

df_comms <- df_comms %>% 
  mutate(prob = ifelse(prob < 1e-7, NA, prob)) %>%
  mutate(prob = ifelse(
    is.na(prob),
    NA,
    (prob - min(prob, na.rm = TRUE)) /
      (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
  ))

df_border_comms <- df_comms %>%
  mutate(has_value = prob > 0.75)  

# Build FDivergence text block
fe_text <- paste0(
  "Evenness = ",
  round(RED_comm_d2$communities$FEvenness, 2),
  collapse = "\n"
)

fr_text <- paste0(
  "Richness = ",
  round(RED_comm_d2$communities$FRichness, 2),
  collapse = "\n"
)

fd_text <- paste0(
  "Divergence = ",
  round(RED_comm_d2$communities$FDivergence, 2),
  collapse = "\n"
)

# Build Redundancy text block
fred_text <- paste0("Redundancy = ", round(FRed_d2$redundancy, 3))

# Combine both
annot_text <- paste(fr_text, fe_text, fd_text, fred_text, sep = "\n\n\n\n")

# Plot
tpd <- ggplot(df_comms, aes(x = Axis.1, y = Axis.2, fill = prob)) +
  geom_raster() +
  geom_contour(
    data = df_border_comms,
    aes(z = as.numeric(has_value)),
    color = "yellow",
    size = 0.1
  ) +
  scale_fill_gradient(
    low = "lightblue",
    high = "darkblue",
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "Beak Length vs Beak Depth",
    x = "Beak Length",
    y = "Beak Depth"
  ) +
  theme_minimal() +
  #xlim(0, 400) + ylim(0, 50) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8)
  ) +
  # Add FDivergence + Redundancy annotation text
  annotate(
    "text",
    x = min(df_comms$Axis.1, na.rm = TRUE),   # use df_comms
    y = max(df_comms$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))


ggsave(paste0("TPC_plot_new/TPDc_beak_comms0922.pdf"), plot = tpd, width = 8, height = 6, units = "in", bg = "transparent")


# Example for a single community
tpd_data_inthigh <- TPDc_divergence_d2_inthigh$data  # Get the first community's data
grid_inthigh <- tpd_data_inthigh$evaluation_grid
prob_inthigh <- TPDc_divergence_d2_inthigh$TPDc$TPDc$Comm.1

# Build data frame for ggplot

df_comms_inthigh <- data.frame(grid_inthigh, prob=prob_inthigh)

df_comms_inthigh <- df_comms_inthigh %>%
  mutate(prob = ifelse(prob < 1e-7, NA, prob)) %>%
  mutate(prob = ifelse(
    is.na(prob),
    NA,
    (prob - min(prob, na.rm = TRUE)) /
      (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
  ))
#df_inthigh <- df_inthigh %>% drop_na(Axis.1, Axis.2, prob)
# Create a mask for the zero/non-zero areas
df_comms_border_inthigh <- df_comms_inthigh %>%
  mutate(has_value = prob > 0.75)  # TRUE if value > 0

# Build FDivergence text block

fr_text <- paste0(
  "Divergence = ",
  round(RED_comm_d2_inthigh$communities$FRichness, 2),
  collapse = "\n"
)

fe_text <- paste0(
  "Evenness = ",
  round(RED_comm_d2_inthigh$communities$FEvenness, 2),
  collapse = "\n"
)

fd_text <- paste0(
  "Divergence=",
  round(RED_comm_d2_inthigh$communities$FDivergence, 2),
  collapse = "\n"
)

# Build Redundancy text block
fred_text <- paste(
  paste0(
         "Redundancy=", round(FRed_d2_inthigh$redundancy, 3)),
  collapse = "\n"
)

# Combine both
annot_text <- paste(fr_text, fe_text, fd_text, fred_text, sep = "\n\n\n\n")

# Plot
tpd_inthigh <- ggplot(df_comms_inthigh, aes(x = Axis.1, y = Axis.2, fill = prob)) +
  geom_raster() +
  geom_contour(
    data = df_comms_border_inthigh,
    aes(z = as.numeric(has_value)),
    color = "red",
    size = 0.4
  ) +
  scale_fill_gradient(
    low = "#FFD6E0",
    high = "#8B004B",
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "Beak Length vs Beak Depth",
    x = "Beak Length",
    y = "Beak Depth"
  ) +
  #xlim(0, 400) + ylim(0, 50) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8)
  ) +
  annotate(
    "text",
    x = min(df_comms_inthigh$Axis.1, na.rm = TRUE),   # safer to use matching df
    y = max(df_comms_inthigh$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))

ggsave(paste0("TPC_plot_new/TPDc_beak_comms_inthigh0922.pdf"), plot = tpd_inthigh, width = 8, height = 6, units = "in", bg = "transparent")



#--------start TPDc and REND calculations for movement trait space reserve level------------------------------
TPDs_d2 <- TPDs(species = traits_with_abund_2comms$species, traits_df2_comms)
RED_d2 <- REND(TPDs = TPDs_d2)



TPDs_d2_inthigh <- TPDs(species = traits_with_abund_2comms_inthigh$species, traits_df2_comms_inthigh)
RED_d2_inthigh <- REND(TPDs = TPDs_d2_inthigh)


df_two_comms <- bind_rows(
  sf_in_any %>% mutate(region = "Historical"),
  sf_in_any_inthigh %>% mutate(region = "Sea Level Rise")
) %>%
  mutate(region = factor(region, levels = c("Historical", "Sea Level Rise"))) %>%
  select(-longitude, -latitude) %>%
  group_by(region) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  select(any_of(unique(traits_with_abund_2comms_inthigh$species)))
df_two_matrix<- as.matrix(df_two_comms)
rownames(df_two_matrix) <- paste0("Comm.", seq_len(nrow(df_two_matrix)))


TPDc_divergence_d2 <- TPDc(TPDs = TPDs_d2, sampUnit = sampUnit)
RED_comm_d2 <- REND(TPDc = TPDc_divergence_d2)

TPDc_divergence_d2_inthigh <- TPDc(TPDs = TPDs_d2_inthigh, sampUnit = sampUnit_inthigh)
RED_comm_d2_inthigh <- REND(TPDc = TPDc_divergence_d2_inthigh)

dissim_d2 <- dissim(TPDs_d2)
dissim_d2_inthigh <- dissim(TPDs_d2_inthigh)

FRed_d2 <- redundancy(TPDc = TPDc_divergence_d2)
FRed_d2_inthigh <- redundancy(TPDc = TPDc_divergence_d2_inthigh)


# Example for a single community
tpd_data_ <- TPDc_divergence_d2$data  # Get the data
grid <- tpd_data_$evaluation_grid
prob <- TPDc_divergence_d2$TPDc$TPDc$Comm.1

# Build data frame for ggplot
df_comms <- data.frame(grid, prob = prob)

df_comms <- df_comms %>% 
  mutate(prob = ifelse(prob < 1e-7, NA, prob)) %>%
  mutate(prob = ifelse(
    is.na(prob),
    NA,
    (prob - min(prob, na.rm = TRUE)) /
      (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
  ))

df_border_comms <- df_comms %>%
  mutate(has_value = prob > 0.75)  

# Build FDivergence text block
fe_text <- paste0(
  "Evenness = ",
  round(RED_comm_d2$communities$FEvenness, 2),
  collapse = "\n"
)

fr_text <- paste0(
  "Richness = ",
  round(RED_comm_d2$communities$FRichness, 2),
  collapse = "\n"
)

fd_text <- paste0(
  "Divergence = ",
  round(RED_comm_d2$communities$FDivergence, 2),
  collapse = "\n"
)

# Build Redundancy text block
fred_text <- paste0("Redundancy = ", round(FRed_d2$redundancy, 3))

# Combine both
annot_text <- paste(fr_text, fe_text, fd_text, fred_text, sep = "\n\n\n\n")

# Plot
tpd <- ggplot(df_comms, aes(x = Axis.1, y = Axis.2, fill = prob)) +
  geom_raster() +
  geom_contour(
    data = df_border_comms,
    aes(z = as.numeric(has_value)),
    color = "yellow",
    size = 0.1
  ) +
  scale_fill_gradient(
    low = "lightblue",
    high = "darkblue",
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "Tail Length vs Hand Wind Index",
    x = "Tail Length",
    y = "Hand Wind Index"
  ) +
  theme_minimal() +
  #xlim(0, 400) + ylim(0, 50) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8)
  ) +
  # Add FDivergence + Redundancy annotation text
  annotate(
    "text",
    x = min(df_comms$Axis.1, na.rm = TRUE),   # use df_comms
    y = max(df_comms$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))


ggsave(paste0("TPC_plot_new/TPDc_movement_comms0922.pdf"), plot = tpd, width = 8, height = 6, units = "in", bg = "transparent")


# Example for a single community
tpd_data_inthigh <- TPDc_divergence_d2_inthigh$data  # Get the first community's data
grid_inthigh <- tpd_data_inthigh$evaluation_grid
prob_inthigh <- TPDc_divergence_d2_inthigh$TPDc$TPDc$Comm.2

# Build data frame for ggplot

df_comms_inthigh <- data.frame(grid_inthigh, prob=prob_inthigh)

df_comms_inthigh <- df_comms_inthigh %>%
  mutate(prob = ifelse(prob < 1e-7, NA, prob)) %>%
  mutate(prob = ifelse(
    is.na(prob),
    NA,
    (prob - min(prob, na.rm = TRUE)) /
      (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
  ))
#df_inthigh <- df_inthigh %>% drop_na(Axis.1, Axis.2, prob)
# Create a mask for the zero/non-zero areas
df_comms_border_inthigh <- df_comms_inthigh %>%
  mutate(has_value = prob > 0.75)  # TRUE if value > 0

# Build FDivergence text block

fr_text <- paste0(
  "Divergence = ",
  round(RED_comm_d2_inthigh$communities$FRichness, 2),
  collapse = "\n"
)

fe_text <- paste0(
  "Evenness = ",
  round(RED_comm_d2_inthigh$communities$FEvenness, 2),
  collapse = "\n"
)

fd_text <- paste0(
  "Divergence=",
  round(RED_comm_d2_inthigh$communities$FDivergence, 2),
  collapse = "\n"
)

# Build Redundancy text block
fred_text <- paste(
  paste0(
    "Redundancy=", round(FRed_d2_inthigh$redundancy, 3)),
  collapse = "\n"
)

# Combine both
annot_text <- paste(fr_text, fe_text, fd_text, fred_text, sep = "\n\n\n\n")

# Plot
tpd_inthigh <- ggplot(df_comms_inthigh, aes(x = Axis.1, y = Axis.2, fill = prob)) +
  geom_raster() +
  geom_contour(
    data = df_comms_border_inthigh,
    aes(z = as.numeric(has_value)),
    color = "red",
    size = 0.4
  ) +
  scale_fill_gradient(
    low = "#FFD6E0",
    high = "#8B004B",
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "Tail Length vs Hand Wing Index",
    x = "Tail Length",
    y = "Hand Wing Depth"
  ) +
  #xlim(0, 400) + ylim(0, 50) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8)
  ) +
  annotate(
    "text",
    x = min(df_comms_inthigh$Axis.1, na.rm = TRUE),   # safer to use matching df
    y = max(df_comms_inthigh$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))

ggsave(paste0("TPC_plot_new/TPDc_movement_comms_inthigh0922.pdf"), plot = tpd_inthigh, width = 8, height = 6, units = "in", bg = "transparent")


# beak plots
TPDs_d2 <- TPDs(species = traits_with_abund_2comms$species, traits_df3_comms)
RED_d2 <- REND(TPDs = TPDs_d2)



TPDs_d2_inthigh <- TPDs(species = traits_with_abund_2comms_inthigh$species, traits_df3_comms_inthigh)
RED_d2_inthigh <- REND(TPDs = TPDs_d2_inthigh)


TPDc_divergence_d2 <- TPDc(TPDs = TPDs_d2, sampUnit = sampUnit)
RED_comm_d2 <- REND(TPDc = TPDc_divergence_d2)

TPDc_divergence_d2_inthigh <- TPDc(TPDs = TPDs_d2_inthigh, sampUnit = sampUnit_inthigh)
RED_comm_d2_inthigh <- REND(TPDc = TPDc_divergence_d2_inthigh)

dissim_d2 <- dissim(TPDs_d2)
dissim_d2_inthigh <- dissim(TPDs_d2_inthigh)

FRed_d2 <- redundancy(TPDc = TPDc_divergence_d2)
FRed_d2_inthigh <- redundancy(TPDc = TPDc_divergence_d2_inthigh)


# Example for a single community
tpd_data_ <- TPDc_divergence_d2$data  # Get the data
grid <- tpd_data_$evaluation_grid
prob <- TPDc_divergence_d2$TPDc$TPDc$Comm.1

# Build data frame for ggplot
df_comms <- data.frame(grid, prob = prob)

df_comms <- df_comms %>% 
  mutate(prob = ifelse(prob < 1e-7, NA, prob)) %>%
  mutate(prob = ifelse(
    is.na(prob),
    NA,
    (prob - min(prob, na.rm = TRUE)) /
      (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
  ))

df_border_comms <- df_comms %>%
  mutate(has_value = prob > 0.75)  

# Build FDivergence text block
fe_text <- paste0(
  "Evenness = ",
  round(RED_comm_d2$communities$FEvenness, 2),
  collapse = "\n"
)

fr_text <- paste0(
  "Richness = ",
  round(RED_comm_d2$communities$FRichness, 2),
  collapse = "\n"
)

fd_text <- paste0(
  "Divergence = ",
  round(RED_comm_d2$communities$FDivergence, 2),
  collapse = "\n"
)

# Build Redundancy text block
fred_text <- paste0("Redundancy = ", round(FRed_d2$redundancy, 3))

# Combine both
annot_text <- paste(fr_text, fe_text, fd_text, fred_text, sep = "\n\n\n\n")

# Plot
tpd <- ggplot(df_comms, aes(x = Axis.1, y = Axis.2, fill = prob)) +
  geom_raster() +
  geom_contour(
    data = df_border_comms,
    aes(z = as.numeric(has_value)),
    color = "yellow",
    size = 0.1
  ) +
  scale_fill_gradient(
    low = "lightblue",
    high = "darkblue",
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "Habitat Density vs Tarsus Length",
    x = "Habitat Density",
    y = "Tarsus Length"
  ) +
  theme_minimal() +
  #xlim(0, 400) + ylim(0, 50) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8)
  ) +
  # Add FDivergence + Redundancy annotation text
  annotate(
    "text",
    x = min(df_comms$Axis.1, na.rm = TRUE),   # use df_comms
    y = max(df_comms$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))


ggsave(paste0("TPC_plot_new/TPDc_size_comms0922.pdf"), plot = tpd, width = 8, height = 6, units = "in", bg = "transparent")


# Example for a single community
tpd_data_inthigh <- TPDc_divergence_d2_inthigh$data  # Get the first community's data
grid_inthigh <- tpd_data_inthigh$evaluation_grid
prob_inthigh <- TPDc_divergence_d2_inthigh$TPDc$TPDc$Comm.2

# Build data frame for ggplot

df_comms_inthigh <- data.frame(grid_inthigh, prob=prob_inthigh)

df_comms_inthigh <- df_comms_inthigh %>%
  mutate(prob = ifelse(prob < 1e-7, NA, prob)) %>%
  mutate(prob = ifelse(
    is.na(prob),
    NA,
    (prob - min(prob, na.rm = TRUE)) /
      (max(prob, na.rm = TRUE) - min(prob, na.rm = TRUE))
  ))
#df_inthigh <- df_inthigh %>% drop_na(Axis.1, Axis.2, prob)
# Create a mask for the zero/non-zero areas
df_comms_border_inthigh <- df_comms_inthigh %>%
  mutate(has_value = prob > 0.75)  # TRUE if value > 0

# Build FDivergence text block

fr_text <- paste0(
  "Divergence = ",
  round(RED_comm_d2_inthigh$communities$FRichness, 2),
  collapse = "\n"
)

fe_text <- paste0(
  "Evenness = ",
  round(RED_comm_d2_inthigh$communities$FEvenness, 2),
  collapse = "\n"
)

fd_text <- paste0(
  "Divergence=",
  round(RED_comm_d2_inthigh$communities$FDivergence, 2),
  collapse = "\n"
)

# Build Redundancy text block
fred_text <- paste(
  paste0(
    "Redundancy=", round(FRed_d2_inthigh$redundancy, 3)),
  collapse = "\n"
)

# Combine both
annot_text <- paste(fr_text, fe_text, fd_text, fred_text, sep = "\n\n\n\n")

# Plot
tpd_inthigh <- ggplot(df_comms_inthigh, aes(x = Axis.1, y = Axis.2, fill = prob)) +
  geom_raster() +
  geom_contour(
    data = df_comms_border_inthigh,
    aes(z = as.numeric(has_value)),
    color = "red",
    size = 0.4
  ) +
  scale_fill_gradient(
    low = "#FFD6E0",
    high = "#8B004B",
    name = "TPD Value",
    na.value = "transparent"
  ) +
  labs(
    title = "Habitat Density vs Tarsus Length",
    x = "Habitat Density",
    y = "Tarsus Length"
  ) +
  #xlim(0, 400) + ylim(0, 50) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8)
  ) +
  annotate(
    "text",
    x = min(df_comms_inthigh$Axis.1, na.rm = TRUE),   # safer to use matching df
    y = max(df_comms_inthigh$Axis.2, na.rm = TRUE),
    label = annot_text,
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))

ggsave(paste0("TPC_plot_new/TPDc_size_comms_inthigh0922.pdf"), plot = tpd_inthigh, width = 8, height = 6, units = "in", bg = "transparent")




