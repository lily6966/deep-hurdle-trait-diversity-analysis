library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(tibble)
library(readr)
setwd("/Users/liyingnceas/GitHub/Deep-hurdle-population-model-exp/lightning_logs")
# --- Inputs ---
# species_df: 7800 rows × 332 species columns (occurrence probabilities)
# env_df:     7800 rows × 28 environmental variables
species_df <- read_csv("inthigh_mean_predictions_month_6.csv")
env_df     <- read_csv("pred_pland_elev_poc_chl_t_trangeintHigh_10year.csv")
label_label_his <- read_csv("historical/label_label_emb_regressor.csv")
label_label_int <- read_csv("intHigh/label_label_emb_regressor.csv")
label_lable_class_his <- read_csv("historical/label_label_emb_classifier.csv")
env_df     <- read_csv("pred_pland_elev_poc_chl_t_trange_10year.csv")
# Keep numeric columns and align rows
species_mat <- species_df %>% select(where(is.numeric)) %>% as.matrix()
env_mat     <- env_df     %>% select(where(is.numeric))  %>% select (-id, -longitude, -eastness_sd, -northness_sd, -eastness_median,
                                                                                    -northness_median, -elevation_median, -latitude, -elevation_sd, -slopeness_sd,
                                                                                    -slopeness_median, -no_value) %>% as.matrix()
stopifnot(nrow(species_mat) == nrow(env_mat))

# Drop zero-variance columns to avoid NA correlations
nzv_species <- apply(species_mat, 2, function(x) sd(x, na.rm = TRUE) > 0)
nzv_env     <- apply(env_mat, 2, function(x) sd(x, na.rm = TRUE) > 0)
species_mat <- species_mat[, nzv_species, drop = FALSE]
env_mat     <- env_mat[, nzv_env, drop = FALSE]

# --- Correlations: env (rows) × species (cols) ---
# Use "spearman" for robustness; switch to "pearson" if preferred
corr <- cor(env_mat, species_mat, method = "spearman", use = "pairwise.complete.obs")
# corr is (#env_vars) × (#species)

# --- Summarize distribution across species per environmental feature ---
corr_summary <- as.data.frame(corr) %>%
  rownames_to_column("env_var") %>%
  pivot_longer(-env_var, names_to = "species", values_to = "rho") %>%
  group_by(env_var) %>%
  summarize(
    mean_r  = mean(rho, na.rm = TRUE),
    median_r = median(rho, na.rm = TRUE),
    sd_r    = sd(rho, na.rm = TRUE),
    iqr_r   = IQR(rho, na.rm = TRUE),
    mean_abs = mean(abs(rho), na.rm = TRUE),
    n       = sum(!is.na(rho)),
    .groups = "drop"
  ) %>%
  mutate(env_var = fct_reorder(env_var, mean_abs, .desc = TRUE))  # order by overall strength

# --- Balloon plot: each point is one env feature ---
dir.create("species_env_corr_outputs", showWarnings = FALSE, recursive = TRUE)

p <- ggplot(corr_summary, aes(x = env_var, y = mean_r)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
  geom_point(aes(size = sd_r, fill = mean_r), shape = 21, color = "grey20", alpha = 0.9) +
  scale_size_continuous(range = c(3, 14), name = "Spread (SD of ρ)") +
  scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#d7301f",
                       midpoint = 0, name = "Mean Spearman ρ") +
  labs(title = "Species–Environment Correlation Summary by Feature",
       x = "Environmental feature (ordered by mean |ρ|)",
       y = "Mean Spearman ρ across species") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

p


ggsave("species_env_corr_outputs/species_env_balloon_summary_historical_9.pdf",
       p, width = 12, height = 7)

# Also save the summary table
write_csv(corr_summary, "species_env_corr_outputs/species_env_correlation_summary.csv")


df_result = read_csv("~/GitHub/gomes-lab-DMVP-DRNets/eBird_entire/results/prb_change_intHigh/proportional_change_with_trait_month9_trophic.csv")


library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# pick your trait column name
tg_col <- if ("TraitGroup" %in% names(df_result)) "TraitGroup" else "traitgroup"

df_long <- df_result %>%
  pivot_longer(
    cols = matches("^\\d+$"),             # site columns "0","1","2",...
    names_to = "site_id",
    values_to = "delta_prob",
    values_drop_na = TRUE
  ) %>%
  filter(is.finite(delta_prob)) %>%
  mutate(
    site_id = as.integer(site_id),
    !!tg_col := fct_reorder(.data[[tg_col]], delta_prob, .fun = median, na.rm = TRUE)
  )


library(ggridges)

levs <- readRDS("~/GitHub/gomes-lab-DMVP-DRNets/eBird_entire/results/prb_change_intHigh/habitat1.rds")
pal  <- readRDS("~/GitHub/gomes-lab-DMVP-DRNets/eBird_entire/results/prb_change_intHigh/habitat1.rds")

pal["Omnivore"] <- "#9C755F"
levs <- levels(df_long[[tg_col]])
# ensure the trait column matches your saved levels
# remove NA trait groups (and any non-finite Δp, optional)
df_long <- df_long %>%
  dplyr::filter(!is.na(.data[[tg_col]])) %>%
  droplevels()

ggplot(df_long, aes(x = delta_prob, y = .data[[tg_col]],
                    height = after_stat(density),
                    fill = .data[[tg_col]])) +
  stat_density_ridges(geom = "density_ridges", scale = 1.2,
                      rel_min_height = 0.001, adjust = 1, color = "grey15") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  theme_minimal(base_size = 12) +
  coord_cartesian(xlim = quantile(df_long$delta_prob, c(0.01, 0.99), na.rm = TRUE)) +
  scale_fill_manual(values = pal, limits = levs, drop = FALSE,
                    name = "Trophic niche", na.translate = FALSE)

  ggsave("ridgelines_traitgroups_month11_trophic0105.pdf",
         width = 9, height = 4, units = "in", dpi = 300,
         bg = "white")  # solid white background for PNG
