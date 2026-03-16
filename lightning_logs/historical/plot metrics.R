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
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# === Step 1: Trait matrix ===
# Each row is a species, each column is a trait
all_metrics <- read_csv("combined_metrics.csv")
count_predictions <-read_csv("multiplied_result.csv") #1-12 pevelant species and 13-24 including rare species
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
traits_metrics_merged <-traits_filtered%>%
  rename(Label = Species2) %>%
  inner_join(all_metrics, by = "Label")

library(ggplot2)

ggplot(traits_metrics_merged, aes(x = Habitat, y = `F1 Score`, fill = Habitat)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Distribution of F1 Scores by Habitat",
    x = "Habitat Category",
    y = "F1 Score"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




r <- ggplot(traits_metrics_merged, aes(x = Habitat, y = `F1 Score`, fill = Habitat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Distribution of F1 Scores by Habitat",
    x = "Habitat Category",
    y = "F1 Score"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # frame
    #plot.background = element_rect(color = "black", fill = "white", linewidth = 1) # outer frame
  )

# 5) Save to PDF
ggsave("F1_by_Habitat.pdf", plot = r, width = 9, height = 6)

count_filtered <- count_predictions
count_filtered[ , 3:ncol(count_filtered)] <- as.data.frame(
  ifelse(
    as.matrix(occurance_binary[ , 3:ncol(occurance_binary)]) == 1,
    as.matrix(count_predictions[ , 3:ncol(count_predictions)]),
    NA
  )
)

count_filtered_summary <- count_filtered %>%
  filter(if_any(-c(longitude, latitude), ~ !is.na(.)))%>%
  mutate(across(everything(), ~replace_na(., 0)))%>%
  summarise(across(everything(), ~sum(.))) %>%
  select(-longitude, -latitude)%>%
  select(where(~ any(. != 0)))

library(dplyr)
library(ggplot2)
library(forcats)

# 1) Species set = intersection of Label values and the column names of count_filtered_summary
species_set <- intersect(traits_metrics_merged$Label, colnames(count_filtered_summary))

# 2) Filter to those species (and clean NAs)
df_filt <- traits_metrics_merged %>%
  filter(Label %in% species_set) %>%
  filter(!is.na(Habitat), !is.na(`F1 Score`))


# 4) Plot (violin + box + jitter), add panel & outer frames
s <- ggplot(df_filt, aes(x = Habitat, y = `F1 Score`, fill = Habitat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Distribution of F1 Scores by Habitat",
    x = "Habitat Category",
    y = "F1 Score"
  ) +
  ylim(0, 1)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # frame
    #plot.background = element_rect(color = "black", fill = "white", linewidth = 1) # outer frame
  )

# 5) Save to PDF
ggsave("F1_by_Habitat_filtered.pdf", plot = s, width = 9, height = 6)



library(dplyr)
library(ggplot2)

# species overlap
species_set <- intersect(traits_metrics_merged$Label, colnames(count_filtered_summary))

# filtered dataframe
df_filt <- traits_metrics_merged %>%
  filter(Label %in% species_set) %>%
  filter(!is.na(Habitat), !is.na(`F1 Score`))

# violin + box to show distribution
p <- ggplot(df_filt, aes(x = Habitat, y = `F1 Score`, fill = Habitat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of F1 Scores by Habitat (filtered species)",
    x = "Habitat Category",
    y = "F1 Score"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.background = element_rect(color = "black", fill = "white", linewidth = 1)
  )

# save to PDF
ggsave("F1_distribution_by_Habitat.pdf", plot = p, width = 9, height = 6)

library(ggplot2)

# Create the plot
p <- ggplot(traits_metrics_merged, aes(x = Habitat, y = `Normalized RMSE (0-1)`, fill = Habitat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Distribution of Normalized RMSE by Habitat",
    x = "Habitat Category",
    y = "Normalized RMSE (0-1)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # frame
    #plot.background = element_rect(color = "black", fill = "white", linewidth = 1) # outer frame
  )

# Save to PDF
ggsave("NRMSE_by_Habitat.pdf", plot = p, width = 8, height = 6)

q <- ggplot(traits_metrics_merged, aes(x = Habitat, y = `RMSE`, fill = Habitat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Distribution of RMSE by Habitat",
    x = "Habitat Category",
    y = "RMSE"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # frame
    #plot.background = element_rect(color = "black", fill = "white", linewidth = 1) # outer frame
  )


# Save to PDF
ggsave("RMSE_by_Habitat.pdf", plot = q, width = 8, height = 6)

library(RcppCNPy)

# Install if needed
install.packages("RcppCNPy")

library(RcppCNPy)

# Direct load of .npy inside .npz won't work, so first unzip
unzip("binary_results.npz", exdir = "npz_extracted")

# Now you'll see one or more .npy files in "npz_extracted/"
# For example, if the array is called "data.npy":
binary_data <- npyLoad("npz_extracted/data.npy")

dim(binary_data)

