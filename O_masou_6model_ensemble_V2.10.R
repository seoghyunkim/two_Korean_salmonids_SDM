####--------- 1. Environment Setup & Package Loading ---------####

# O. masou masou -- six-model ensemble SDM, V2.10

sessionInfo()

rm(list = ls())

# Load packages
library(ENMeval)   # MaxEnt tuning
library(sf)        
library(usdm)      # VIF
library(ecospat)   # MaxSSS threshold
library(dismo)     # MaxEnt fitting
library(caret)     # RF, GBM, and ANN tuning
library(ranger)    # random forest
library(gbm)       # GBM  
library(nnet)      # ANN
library(xgboost)   # xgboost
library(lightgbm)  # light GBM
library(pROC)      # AUC
library(ggpubr)
library(ggspatial)
library(tidyverse)

set.seed(777)

dir.create("./data_outcome/O_masou_6model_ensemble_V2.10", recursive = TRUE, showWarnings = FALSE)



####--------- 2. Data Loading & Preprocessing ---------####

# Occurrence records are already assigned to stream segments.
presencedata <- st_read("./data_raw/O_masou_occ_corrected_clean_snap_join.shp") %>%
  mutate(BBSNCD = as.character(BBSNCD),
         SBSNCD = as.character(SBSNCD),
         str_ID = as.character(str_ID)) %>%
  filter(BBSNCD %in% c("10", "13", "24")) %>% # select target large basin
  st_transform(4326) %>%
  mutate(slp_CV = std_stslp / mean_stslp) %>%
  rename(ws_area = wsarea, elev_mean = meanstdem, slp_mean = mean_stslp,
         st_power_mean = mean_stpo, urban = per_ct_urb,
         agri = per_ct_agr, forest = per_ct_for) %>%
  distinct(str_ID, .keep_all = TRUE) %>%
  filter(!is.na(p_bio01))


# Large basins are selected from the retained occurrence records.
target_large_basins <- unique(na.omit(presencedata$BBSNCD))


# Stream network
streamnetwork <- st_read("./data_raw/epsg4326_catchment_ver1.2_mkp2.shp") %>%
  mutate(BBSNCD = as.character(BBSNCD),
         SBSNCD = as.character(SBSNCD),
         str_ID = as.character(str_ID)) %>%
  filter(BBSNCD %in% c("10", "13", "24")) %>% # select target large basin
  filter(BBSNCD %in% target_large_basins) %>% # select large basins with occurrence records
  st_transform(4326) %>%
  mutate(slp_CV = std_stslp / mean_stslp) %>% # calculate segment CV
  rename(ws_area = wsarea, elev_mean = meanstdem, slp_mean = mean_stslp,
         st_power_mean = mean_stpo, urban = per_ct_urb,
         agri = per_ct_agr, forest = per_ct_for) %>%
  distinct(str_ID, .keep_all = TRUE) %>%
  filter(!is.na(p_bio01))


# Accessible area: standard basins containing occurrence records
accessible_basins <- unique(presencedata$SBSNCD)

streamnetwork2 <- streamnetwork %>%
  filter(SBSNCD %in% accessible_basins)

streamnetwork3 <- streamnetwork2 %>%
  st_drop_geometry()

bgdata <- streamnetwork2 %>%
  filter(!str_ID %in% presencedata$str_ID)


# Load large basin boundaries
korea_watershed <- st_read("./data_raw/WKMBBSN.shp") %>%
  mutate(BBSNCD = as.character(BBSNCD)) %>%
  st_transform(4326)

map_watershed <- korea_watershed %>%
  filter(BBSNCD %in% target_large_basins)

map_streamnetwork <- streamnetwork %>%
  filter(BBSNCD %in% target_large_basins)



####--------- 3-1. Predictor Variables & VIF ---------####

# Seven predictors used for both species
# BIO05 = Maximum temperature of the warmest month
# BIO13 = Precipitation of the wettest month
env_vars <- c("slp_mean", "slp_CV", "ws_area", "forest", "weir_dens",
              "p_bio05", "p_bio13")

vif_input <- streamnetwork3 %>%
  dplyr::select(all_of(env_vars)) %>%
  mutate(across(everything(), as.numeric)) %>%
  filter(if_all(everything(), is.finite)) %>%
  as.data.frame()

vif_result <- as.data.frame(usdm::vif(vif_input))

vif_result

write.csv(
  vif_result,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_VIF_result.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 3-2. Pearson Correlation & Significance Test ---------####

library(GGally)  # used for correlation matrix visualization

# Select variables for correlation analysis
cor1 <- streamnetwork3 %>%
  dplyr::select(elev_mean, slp_mean, slp_CV, ws_area,
                forest, agri, weir_dens,
                p_bio01, p_bio02, p_bio03, p_bio04, p_bio05,
                p_bio06, p_bio07, p_bio08, p_bio09, p_bio10,
                p_bio11, p_bio12, p_bio13, p_bio14, p_bio15,
                p_bio16, p_bio17, p_bio18, p_bio19) %>%
  rename_with(
    ~ gsub("^p_bio", "Bio", ., ignore.case = TRUE),
    starts_with("p_bio")
  )


# Variable order
vars <- c("elev_mean", "ws_area", "slp_mean", "slp_CV",
          "forest", "agri", "weir_dens",
          "Bio01", "Bio02", "Bio03", "Bio04", "Bio05",
          "Bio06", "Bio07", "Bio08", "Bio09", "Bio10",
          "Bio11", "Bio12", "Bio13", "Bio14", "Bio15",
          "Bio16", "Bio17", "Bio18", "Bio19")

x <- cor1 %>%
  dplyr::select(all_of(vars))

# Calculate correlation coefficients and p-values
rmat <- cor(x, use = "complete.obs")
pmat <- sapply(seq_along(vars), function(i) {
  sapply(seq_along(vars), function(j) {
    cor.test(x[[i]], x[[j]])$p.value
  })
})

df <- expand.grid(Var1 = vars, Var2 = vars) %>%
  mutate(
    i = match(Var1, vars),
    j = match(Var2, vars),
    r = rmat[cbind(i, j)],
    p = pmat[cbind(i, j)],
    sig = case_when(
      p < .001 ~ "***",
      p < .01  ~ "**",
      p < .05  ~ "*",
      TRUE     ~ ""
    ),
    Var1 = factor(Var1, levels = vars),
    Var2 = factor(Var2, levels = vars)
  ) %>%
  filter(i <= j)

cor_fig <- ggplot(df, aes(Var1, Var2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4, vjust = 0.6) +
  geom_text(aes(label = sig), size = 4, vjust = -0.6) +

  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +

  coord_fixed() +

  scale_x_discrete(
    position = "bottom",
    expand = c(0, 0),
    drop = FALSE
  ) +

  scale_y_discrete(
    limits = rev(vars),
    expand = c(0, 0),
    drop = FALSE
  ) +

  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    panel.grid = element_blank(),
    plot.caption = element_text(
      size = 10,
      hjust = 0,
      margin = margin(t = 8),
      colour = "black"
    )
  ) +

  labs(caption = "Significance levels: * p < 0.05, ** p < 0.01, *** p < 0.001")

cor_fig

# Save figures
ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_Pearson_correlation.png",
  plot = cor_fig,
  width = 16,
  height = 15,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_Pearson_correlation.svg",
  plot = cor_fig,
  width = 16,
  height = 15,
  device = grDevices::svg,
  bg = "white"
)



####--------- 4. Occurrence & Background Data Preparation ---------####

# One complete record per occurrence segment
occ <- presencedata %>%
  st_drop_geometry() %>%
  dplyr::select(str_ID, SBSNCD, longitude = mean_lon, latitude = mean_lat,
                all_of(env_vars)) %>%
  mutate(across(all_of(env_vars), as.numeric)) %>%
  na.omit() %>%
  filter(if_all(all_of(env_vars), is.finite)) %>%
  distinct(str_ID, .keep_all = TRUE) %>%
  as_tibble()


# Candidate background segments
bgdata2 <- bgdata %>%
  st_drop_geometry() %>%
  dplyr::select(str_ID, SBSNCD, longitude = mean_lon, latitude = mean_lat,
                all_of(env_vars)) %>%
  mutate(across(all_of(env_vars), as.numeric)) %>%
  na.omit() %>%
  filter(if_all(all_of(env_vars), is.finite)) %>%
  distinct(str_ID, .keep_all = TRUE) %>%
  as_tibble()


# MaxEnt background: up to 10 segments per presence
samplesize <- min(nrow(bgdata2), nrow(occ) * 10)

set.seed(777)
bg <- bgdata2[sample(seq_len(nrow(bgdata2)), samplesize), ]

write.csv(
  occ,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_presence_segments.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  bg,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_background_segments.csv",
  row.names = FALSE,
  na = "NA"
)


sample_counts <- data.frame(
  n_presence = nrow(occ),
  n_presence_basins = n_distinct(occ$SBSNCD),
  n_available_background = nrow(bgdata2),
  n_sampled_background = nrow(bg)
)

sample_counts

write.csv(
  sample_counts,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_sample_counts.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 5-1. Standard-basin Spatial Cross-validation ---------####

# Assign all records from the same SBSNCD to one fold.
n_folds <- 5

pa_model_data <- bind_rows(
  occ %>% mutate(presence = factor("X1", levels = c("X0", "X1")),
                 source = "presence"),
  bg %>% mutate(presence = factor("X0", levels = c("X0", "X1")),
                source = "background")
) %>%
  distinct(str_ID, .keep_all = TRUE) %>%
  mutate(spatial_group_id = paste0("SBSNCD_", SBSNCD))

basin_summary <- pa_model_data %>%
  group_by(spatial_group_id) %>%
  summarise(
    n_presence = sum(presence == "X1"),
    n_background = sum(presence == "X0"),
    n_total = n(),
    .groups = "drop"
  )


# Balance presence counts first, then background counts.
best_score <- Inf

for (trial in seq_len(5000)) {
  set.seed(777 + trial)

  basin_order <- basin_summary %>%
    mutate(random_order = runif(n())) %>%
    arrange(desc(n_presence > 0), desc(n_presence), random_order)

  fold_counts <- data.frame(
    fold_id = seq_len(n_folds),
    n_presence = 0,
    n_background = 0,
    n_total = 0
  )

  assignment <- integer(nrow(basin_order))

  for (i in seq_len(nrow(basin_order))) {
    basin_i <- basin_order[i, ]

    if (basin_i$n_presence > 0) {
      candidates <- fold_counts %>%
        filter(n_presence == min(n_presence)) %>%
        arrange(n_total)
    } else {
      candidates <- fold_counts %>%
        filter(n_background == min(n_background)) %>%
        arrange(n_total)
    }

    candidate_folds <- candidates %>%
      filter(n_total == min(n_total)) %>%
      pull(fold_id)

    chosen_fold <- candidate_folds[sample.int(length(candidate_folds), 1)]
    assignment[i] <- chosen_fold

    fold_counts$n_presence[chosen_fold] <-
      fold_counts$n_presence[chosen_fold] + basin_i$n_presence
    fold_counts$n_background[chosen_fold] <-
      fold_counts$n_background[chosen_fold] + basin_i$n_background
    fold_counts$n_total[chosen_fold] <-
      fold_counts$n_total[chosen_fold] + basin_i$n_total
  }

  if (all(fold_counts$n_presence > 0) && all(fold_counts$n_background > 0)) {
    fold_score <-
      2.0 * sd(fold_counts$n_presence) / (mean(fold_counts$n_presence) + 1) +
      1.0 * sd(fold_counts$n_background) / (mean(fold_counts$n_background) + 1) +
      0.5 * sd(fold_counts$n_total) / (mean(fold_counts$n_total) + 1)

    if (fold_score < best_score) {
      best_score <- fold_score
      fold_assignment <- data.frame(
        spatial_group_id = basin_order$spatial_group_id,
        fold_id = assignment
      )
    }
  }
}

pa_model_data <- pa_model_data %>%
  left_join(fold_assignment, by = "spatial_group_id")

fold_summary <- pa_model_data %>%
  group_by(fold_id) %>%
  summarise(
    n_presence = sum(presence == "X1"),
    n_background = sum(presence == "X0"),
    n_basins = n_distinct(SBSNCD),
    .groups = "drop"
  )

fold_summary

write.csv(
  fold_summary,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_spatial_CV_fold_summary.csv",
  row.names = FALSE,
  na = "NA"
)

fold_records <- pa_model_data %>%
  dplyr::select(str_ID, SBSNCD, source, presence, longitude, latitude,
                spatial_group_id, fold_id)

write.csv(
  fold_records,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_spatial_CV_fold_assignment.csv",
  row.names = FALSE,
  na = "NA"
)


# Separate occurrence and background records for model fitting.
occ_folded <- pa_model_data %>%
  filter(presence == "X1") %>%
  dplyr::select(-presence, -source)

bg_folded <- pa_model_data %>%
  filter(presence == "X0") %>%
  dplyr::select(-presence, -source)



####--------- 5-2. Occurrence & Spatial Fold Maps ---------####

# Map extent
map_bbox <- st_bbox(map_watershed)
map_xpad <- as.numeric(map_bbox["xmax"] - map_bbox["xmin"]) * 0.03
map_ypad <- as.numeric(map_bbox["ymax"] - map_bbox["ymin"]) * 0.03

coord_map <- coord_sf(
  xlim = c(map_bbox["xmin"] - map_xpad, map_bbox["xmax"] + map_xpad),
  ylim = c(map_bbox["ymin"] - map_ypad, map_bbox["ymax"] + map_ypad),
  expand = FALSE
)

point_sf <- st_as_sf(
  pa_model_data,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)


# Occurrence map
occurrence_plot <- ggplot() +
  geom_sf(data = map_watershed, fill = NA, colour = "gray30", linewidth = 0.4) +
  geom_sf(data = map_streamnetwork, colour = "gray80", linewidth = 0.25) +
  geom_sf(data = point_sf[point_sf$presence == "X1", ],
          colour = "#d7191c", size = 1.5) +
  coord_map +
  theme_void(base_size = 11)

occurrence_plot

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_occurrence_plot.png",
  plot = occurrence_plot,
  width = 7,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_occurrence_plot.svg",
  plot = occurrence_plot,
  width = 7,
  height = 6,
  device = grDevices::svg,
  bg = "white"
)


# Point map: draw background first and presence on top.
bg_data_fig <- point_sf %>%
  filter(source == "background")

occ_data_fig <- point_sf %>%
  filter(source == "presence")

occ_bg_plot <- ggplot() +
  geom_sf(data = map_watershed, fill = NA, colour = "gray30", linewidth = 0.4) +
  geom_sf(data = map_streamnetwork, colour = "gray80", linewidth = 0.25) +
  geom_sf(data = bg_data_fig, aes(colour = source), size = 1.1) +
  geom_sf(data = occ_data_fig, aes(colour = source), size = 1.5) +
  scale_colour_manual(
    values = c(background = "blue", presence = "red"),
    breaks = c("presence", "background"),
    labels = c("Presence", "Background"),
    name = "Records"
  ) +
  coord_map +
  ggspatial::annotation_scale(
    location = "bl",
    width_hint = 0.22,
    unit_category = "metric",
    text_cex = 0.7,
    line_width = 0.35,
    pad_x = grid::unit(0.18, "cm"),
    pad_y = grid::unit(0.18, "cm")
  ) +
  theme_void(base_size = 11)

occ_bg_plot

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_occurrence_background_plot.png",
  plot = occ_bg_plot,
  width = 7,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_occurrence_background_plot.svg",
  plot = occ_bg_plot,
  width = 7,
  height = 6,
  device = grDevices::svg,
  bg = "white"
)


# Segment map: use the same sampled records and their original channel geometry.
# Blue segments are sampled background, not confirmed absences or unsuitable habitat.
bg_segment_fig <- map_streamnetwork %>%
  filter(str_ID %in% bg$str_ID) %>%
  mutate(source = "background")

occ_segment_fig <- map_streamnetwork %>%
  filter(str_ID %in% occ$str_ID) %>%
  mutate(source = "presence")

occ_bg_segment_plot <- ggplot() +
  geom_sf(data = map_watershed, fill = NA, colour = "gray30", linewidth = 0.4) +
  geom_sf(data = map_streamnetwork, colour = "gray80", linewidth = 0.25) +
  geom_sf(data = bg_segment_fig, aes(colour = source), linewidth = 0.7) +
  geom_sf(data = occ_segment_fig, aes(colour = source), linewidth = 1.0) +
  scale_colour_manual(
    values = c(background = "blue", presence = "red"),
    breaks = c("presence", "background"),
    labels = c("Presence", "Background"),
    name = "Records"
  ) +
  coord_map +
  ggspatial::annotation_scale(
    location = "bl",
    width_hint = 0.22,
    unit_category = "metric",
    text_cex = 0.7,
    line_width = 0.35,
    pad_x = grid::unit(0.18, "cm"),
    pad_y = grid::unit(0.18, "cm")
  ) +
  theme_void(base_size = 11)

occ_bg_segment_plot

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_occurrence_background_segment_plot.png",
  plot = occ_bg_segment_plot,
  width = 7,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_occurrence_background_segment_plot.svg",
  plot = occ_bg_segment_plot,
  width = 7,
  height = 6,
  device = grDevices::svg,
  bg = "white"
)


# Spatial fold map
fold_segments <- pa_model_data %>%
  st_drop_geometry() %>%
  dplyr::select(str_ID, fold_id) %>%
  distinct()

fold_map_data <- map_streamnetwork %>%
  left_join(fold_segments, by = "str_ID")

fold_plot <- ggplot() +
  geom_sf(data = map_streamnetwork, colour = "gray88", linewidth = 0.22) +
  geom_sf(data = fold_map_data %>% filter(!is.na(fold_id)),
          aes(colour = factor(fold_id)), linewidth = 0.45) +
  geom_sf(data = map_watershed, fill = NA, colour = "gray30", linewidth = 0.4) +
  coord_map +
  ggspatial::annotation_scale(
    location = "bl",
    width_hint = 0.22,
    text_cex = 0.7,
    line_width = 0.35,
    pad_x = grid::unit(0.18, "cm"),
    pad_y = grid::unit(0.18, "cm")
  ) +
  scale_colour_brewer(palette = "Set1", name = "Spatial CV fold") +
  theme_void(base_size = 11) +
  theme(
    legend.position = "right"
  )

fold_plot

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_spatial_CV_fold_map.png",
  plot = fold_plot,
  width = 7.5,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_spatial_CV_fold_map.svg",
  plot = fold_plot,
  width = 7.5,
  height = 6,
  device = grDevices::svg,
  bg = "white"
)



####--------- 6-1. MaxEnt Model Tuning ---------####

# Coordinates followed by predictor variables
occ_maxent <- occ_folded %>%
  dplyr::select(longitude, latitude, all_of(env_vars)) %>%
  as.data.frame()

bg_maxent <- bg_folded %>%
  dplyr::select(longitude, latitude, all_of(env_vars)) %>%
  as.data.frame()

maxent_user_groups <- list(
  occs.grp = occ_folded$fold_id,
  bg.grp = bg_folded$fold_id
)

# Feature classes and regularization multipliers
tune.args <- list(
  fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
  rm = c(1, 2, 3, 4, 5, 7, 10)
)

mod_res <- ENMeval::ENMevaluate(
  occs = occ_maxent,
  bg = bg_maxent,
  algorithm = "maxent.jar",
  tune.args = tune.args,
  partitions = "user",
  user.grp = maxent_user_groups,
  other.settings = list(pred.type = "cloglog", validation.bg = "partition"),
  parallel = FALSE
)

maxent_tuning <- ENMeval::eval.results(mod_res)

write.csv(
  maxent_tuning,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_MaxEnt_tuning_results.csv",
  row.names = FALSE,
  na = "NA"
)


# Select the lowest AICc.
sel_opt <- maxent_tuning %>%
  filter(is.finite(AICc)) %>%
  arrange(AICc) %>%
  slice(1)

selected_fc <- as.character(sel_opt$fc[1])
selected_rm <- as.numeric(sel_opt$rm[1])

sel_opt

write.csv(
  sel_opt,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_MaxEnt_selected_parameters.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 6-2. MaxEnt Spatial Evaluation & Final Model ---------####

# Apply the selected feature classes and regularization multiplier.
maxent_args <- c(
  paste0("betamultiplier=", selected_rm),
  "autofeature=false",
  paste0("linear=", tolower(grepl("L", selected_fc))),
  paste0("quadratic=", tolower(grepl("Q", selected_fc))),
  paste0("hinge=", tolower(grepl("H", selected_fc))),
  paste0("product=", tolower(grepl("P", selected_fc))),
  paste0("threshold=", tolower(grepl("T", selected_fc))),
  "outputformat=cloglog",
  "jackknife=false",
  "responsecurves=false"
)


# Predict each validation fold using the other four folds.
maxent_oof_list <- vector("list", n_folds)

for (f in seq_len(n_folds)) {
  train_set <- pa_model_data %>% filter(fold_id != f)
  valid_set <- pa_model_data %>% filter(fold_id == f)

  fold_path <- paste0("./data_outcome/O_masou_6model_ensemble_V2.10/MaxEnt_fold_", f)
  dir.create(fold_path, recursive = TRUE, showWarnings = FALSE)

  model_f <- dismo::maxent(
    x = as.data.frame(train_set[, env_vars]),
    p = as.integer(train_set$presence == "X1"),
    args = maxent_args,
    path = fold_path
  )

  pred_f <- predict(
    model_f,
    x = as.data.frame(valid_set[, env_vars]),
    args = "outputformat=cloglog"
  )

  maxent_oof_list[[f]] <- valid_set %>%
    transmute(str_ID, fold_id,
              obs = as.integer(presence == "X1"),
              prob_maxent = as.numeric(pred_f))
}

maxent_oof_pred <- bind_rows(maxent_oof_list)

write.csv(
  maxent_oof_pred,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_MaxEnt_spatial_OOF_predictions.csv",
  row.names = FALSE,
  na = "NA"
)


# Fit the selected model using all occurrence and background records.
dir.create("./data_outcome/O_masou_6model_ensemble_V2.10/MaxEnt_final", recursive = TRUE, showWarnings = FALSE)

mod_maxent <- dismo::maxent(
  x = as.data.frame(pa_model_data[, env_vars]),
  p = as.integer(pa_model_data$presence == "X1"),
  args = maxent_args,
  path = "./data_outcome/O_masou_6model_ensemble_V2.10/MaxEnt_final"
)

saveRDS(mod_maxent, "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_MaxEnt_final.rds")



####--------- 7. Repeated Pseudoabsence Datasets ---------####

# Ten datasets; the first five are used for hyperparameter tuning.
n_pa_reps <- 10
n_tuning_reps <- 5

ml_pa_dataset_list <- vector("list", n_pa_reps)

set.seed(777)
for (r in seq_len(n_pa_reps)) {
  sampled_bg <- list()

  for (f in sort(unique(occ_folded$fold_id))) {
    occ_f <- occ_folded %>% filter(fold_id == f)
    bg_f <- bg_folded %>% filter(fold_id == f)

    # Up to three pseudoabsences per presence within each fold
    n_bg <- min(nrow(bg_f), nrow(occ_f) * 3)
    sampled_bg[[as.character(f)]] <- slice_sample(bg_f, n = n_bg, replace = FALSE)
  }

  ml_pa_dataset_list[[r]] <- bind_rows(
    occ_folded %>% mutate(presence = factor("X1", levels = c("X0", "X1")),
                         source = "presence"),
    bind_rows(sampled_bg) %>%
      mutate(presence = factor("X0", levels = c("X0", "X1")),
             source = "pseudoabsence")
  ) %>%
    mutate(pa_rep = r, row_key = row_number())
}

pa_records <- bind_rows(ml_pa_dataset_list) %>%
  dplyr::select(pa_rep, str_ID, SBSNCD, fold_id, presence)

pa_summary <- pa_records %>%
  count(pa_rep, presence)

pa_summary

write.csv(
  pa_records,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_repeated_PA_records.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  pa_summary,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_repeated_PA_summary.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 8. Random Forest Hyperparameter Tuning ---------####

# AUC and binary metrics; a supplied threshold is kept fixed.
calculate_metrics <- function(prob, obs, threshold = NULL) {
  metric_input <- data.frame(prob = as.numeric(prob), obs = as.integer(obs)) %>%
    filter(is.finite(prob), !is.na(obs))

  if (is.null(threshold)) {
    threshold <- as.numeric(
      ecospat::ecospat.max.tss(metric_input$prob, metric_input$obs)$max.threshold
    )
  }

  pred <- as.integer(metric_input$prob >= threshold)

  TP <- sum(pred == 1 & metric_input$obs == 1)
  TN <- sum(pred == 0 & metric_input$obs == 0)
  FP <- sum(pred == 1 & metric_input$obs == 0)
  FN <- sum(pred == 0 & metric_input$obs == 1)

  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)

  roc_result <- pROC::roc(
    response = metric_input$obs,
    predictor = metric_input$prob,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  data.frame(
    maxsss_thr = threshold,
    AUC = as.numeric(pROC::auc(roc_result)),
    TSS = sensitivity + specificity - 1,
    sensitivity = sensitivity,
    specificity = specificity,
    TP = TP, TN = TN, FP = FP, FN = FN
  )
}


# caret tuning criterion: TSS, with AUC retained for model selection.
tss_summary <- function(data, lev = NULL, model = NULL) {
  metrics <- calculate_metrics(data$X1, as.integer(data$obs == "X1"))

  result <- c(
    TSS = metrics$TSS,
    AUC = metrics$AUC,
    Sensitivity = metrics$sensitivity,
    Specificity = metrics$specificity
  )

  result[is.na(result)] <- 0
  return(result)
}


foreach::registerDoSEQ()

rf_grid <- expand.grid(
  mtry = sort(unique(pmax(1, pmin(length(env_vars), c(2, 3, 4, 5, 6))))),
  splitrule = c("gini", "hellinger"),
  min.node.size = c(5, 10, 15)
)

extra_grid <- expand.grid(
  num_trees = c(500, 1000),
  sample.fraction = c(0.632, 0.70)
)

rf_tuning_all <- list()

for (r in seq_len(n_tuning_reps)) {
  message("RF tuning with repeated pseudoabsence dataset: ", r)

  train_set_r <- ml_pa_dataset_list[[r]] %>%
    dplyr::select(presence, all_of(env_vars), fold_id)

  index <- list()
  indexOut <- list()

  for (f in seq_len(n_folds)) {
    index[[paste0("Fold", f)]] <- which(train_set_r$fold_id != f)
    indexOut[[paste0("Fold", f)]] <- which(train_set_r$fold_id == f)
  }

  control_ml_r <- caret::trainControl(
    method = "cv",
    number = n_folds,
    index = index,
    indexOut = indexOut,
    verboseIter = FALSE,
    classProbs = TRUE,
    summaryFunction = tss_summary,
    savePredictions = "final",
    sampling = NULL,
    allowParallel = FALSE
  )

  train_set_ml_r <- train_set_r %>%
    dplyr::select(presence, all_of(env_vars))

  for (i in seq_len(nrow(extra_grid))) {
    message("RF repetition ", r, ": grid ", i, "/", nrow(extra_grid))

    rf_model_r <- caret::train(
      presence ~ .,
      data = train_set_ml_r,
      method = "ranger",
      num.trees = extra_grid$num_trees[i],
      sample.fraction = extra_grid$sample.fraction[i],
      importance = "permutation",
      trControl = control_ml_r,
      metric = "TSS",
      maximize = TRUE,
      tuneGrid = rf_grid
    )

    rf_tuning_all[[length(rf_tuning_all) + 1]] <- rf_model_r$results %>%
      mutate(
        pa_rep = r,
        num_trees = extra_grid$num_trees[i],
        sample.fraction = extra_grid$sample.fraction[i]
      )
  }
}

rf_tuning_result <- bind_rows(rf_tuning_all)

rf_tuning_summary <- rf_tuning_result %>%
  group_by(mtry, splitrule, min.node.size, num_trees, sample.fraction) %>%
  summarise(
    TSS_mean = mean(TSS, na.rm = TRUE),
    TSS_sd = sd(TSS, na.rm = TRUE),
    AUC_mean = mean(AUC, na.rm = TRUE),
    Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
    Specificity_mean = mean(Specificity, na.rm = TRUE),
    n_tuning_reps = n_distinct(pa_rep),
    .groups = "drop"
  ) %>%
  arrange(desc(TSS_mean), desc(AUC_mean))

best_rf_params <- rf_tuning_summary %>%
  slice(1)

best_rf_params

write.csv(
  rf_tuning_result,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_RF_tuning_raw.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  rf_tuning_summary,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_RF_tuning_summary.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  best_rf_params,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_RF_selected_parameters.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 9. GBM Hyperparameter Tuning ---------####

gbmGrid <- expand.grid(
  n.trees = c(500, 1000, 1500, 2000),
  interaction.depth = c(2, 3, 5),
  shrinkage = c(0.05, 0.01, 0.005),
  n.minobsinnode = c(10, 15, 20)
)

gbm_tuning_all <- list()


for (r in seq_len(n_tuning_reps)) {
  message("GBM tuning with repeated pseudoabsence dataset: ", r)

  train_set_r <- ml_pa_dataset_list[[r]] %>%
    dplyr::select(presence, all_of(env_vars), fold_id)

  index <- list()
  indexOut <- list()

  for (f in seq_len(n_folds)) {
    index[[paste0("Fold", f)]] <- which(train_set_r$fold_id != f)
    indexOut[[paste0("Fold", f)]] <- which(train_set_r$fold_id == f)
  }

  control_ml_r <- caret::trainControl(
    method = "cv",
    number = n_folds,
    index = index,
    indexOut = indexOut,
    verboseIter = FALSE,
    classProbs = TRUE,
    summaryFunction = tss_summary,
    savePredictions = "final",
    sampling = NULL,
    allowParallel = FALSE
  )

  train_set_ml_r <- train_set_r %>%
    dplyr::select(presence, all_of(env_vars))

  gbm_model_r <- caret::train(
    presence ~ .,
    data = train_set_ml_r,
    method = "gbm",
    trControl = control_ml_r,
    verbose = FALSE,
    metric = "TSS",
    maximize = TRUE,
    tuneGrid = gbmGrid,
    bag.fraction = 0.70
  )

  gbm_tuning_all[[r]] <- gbm_model_r$results %>%
    mutate(pa_rep = r)
}

gbm_tuning_result <- bind_rows(gbm_tuning_all)

gbm_tuning_summary <- gbm_tuning_result %>%
  group_by(n.trees, interaction.depth, shrinkage, n.minobsinnode) %>%
  summarise(
    TSS_mean = mean(TSS, na.rm = TRUE),
    TSS_sd = sd(TSS, na.rm = TRUE),
    AUC_mean = mean(AUC, na.rm = TRUE),
    Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
    Specificity_mean = mean(Specificity, na.rm = TRUE),
    n_tuning_reps = n_distinct(pa_rep),
    .groups = "drop"
  ) %>%
  arrange(desc(TSS_mean), desc(AUC_mean))

best_gbm_params <- gbm_tuning_summary %>%
  slice(1)

best_gbm_params

write.csv(
  gbm_tuning_result,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_GBM_tuning_raw.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  gbm_tuning_summary,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_GBM_tuning_summary.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  best_gbm_params,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_GBM_selected_parameters.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 10. ANN Hyperparameter Tuning ---------####

ann_grid <- expand.grid(
  size = c(2, 4, 6, 8),
  decay = c(0.001, 0.01, 0.1)
)

ann_tuning_all <- list()


for (r in seq_len(n_tuning_reps)) {
  message("ANN tuning with repeated pseudoabsence dataset: ", r)

  train_set_r <- ml_pa_dataset_list[[r]] %>%
    dplyr::select(presence, all_of(env_vars), fold_id)

  index <- list()
  indexOut <- list()

  for (f in seq_len(n_folds)) {
    index[[paste0("Fold", f)]] <- which(train_set_r$fold_id != f)
    indexOut[[paste0("Fold", f)]] <- which(train_set_r$fold_id == f)
  }

  control_ml_r <- caret::trainControl(
    method = "cv",
    number = n_folds,
    index = index,
    indexOut = indexOut,
    verboseIter = FALSE,
    classProbs = TRUE,
    summaryFunction = tss_summary,
    savePredictions = "final",
    sampling = NULL,
    allowParallel = FALSE
  )

  train_set_ml_r <- train_set_r %>%
    dplyr::select(presence, all_of(env_vars))

  ann_model_r <- caret::train(
    presence ~ .,
    data = train_set_ml_r,
    method = "nnet",
    trControl = control_ml_r,
    metric = "TSS",
    maximize = TRUE,
    tuneGrid = ann_grid,
    preProcess = c("center", "scale"),
    trace = FALSE,
    maxit = 500,
    MaxNWts = 10000
  )

  ann_tuning_all[[r]] <- ann_model_r$results %>%
    mutate(pa_rep = r)
}

ann_tuning_result <- bind_rows(ann_tuning_all)

ann_tuning_summary <- ann_tuning_result %>%
  group_by(size, decay) %>%
  summarise(
    TSS_mean = mean(TSS, na.rm = TRUE),
    TSS_sd = sd(TSS, na.rm = TRUE),
    AUC_mean = mean(AUC, na.rm = TRUE),
    Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
    Specificity_mean = mean(Specificity, na.rm = TRUE),
    n_tuning_reps = n_distinct(pa_rep),
    .groups = "drop"
  ) %>%
  arrange(desc(TSS_mean), desc(AUC_mean))

best_ann_params <- ann_tuning_summary %>%
  slice(1)

best_ann_params

write.csv(
  ann_tuning_result,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ANN_tuning_raw.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  ann_tuning_summary,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ANN_tuning_summary.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  best_ann_params,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ANN_selected_parameters.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 11. XGBoost Hyperparameter Tuning ---------####

xgb_grid <- expand.grid(
  nrounds = c(100, 300),
  max_depth = c(2, 3, 5),
  eta = c(0.01, 0.05),
  gamma = c(0, 1),
  colsample_bytree = c(0.8),
  min_child_weight = c(1, 5),
  subsample = c(0.7)
)

xgb_tuning_all <- list()

for (r in seq_len(n_tuning_reps)) {
  message("XGBoost tuning with repeated pseudoabsence dataset: ", r)

  train_set_r <- ml_pa_dataset_list[[r]] %>%
    dplyr::select(presence, all_of(env_vars), fold_id)

  fold_values_xgb <- sort(unique(train_set_r$fold_id))

  for (g in seq_len(nrow(xgb_grid))) {
    message("  XGBoost grid ", g, " / ", nrow(xgb_grid))

    grid_i <- xgb_grid[g, ]
    xgb_fold_pred_list <- list()

    for (f in fold_values_xgb) {
      train_fold_i <- train_set_r %>% filter(fold_id != f)
      valid_fold_i <- train_set_r %>% filter(fold_id == f)

      train_y <- ifelse(train_fold_i$presence == "X1", 1, 0)
      valid_y <- ifelse(valid_fold_i$presence == "X1", 1, 0)

      train_x <- as.matrix(as.data.frame(train_fold_i[, env_vars]))
      valid_x <- as.matrix(as.data.frame(valid_fold_i[, env_vars]))

      dtrain <- xgboost::xgb.DMatrix(data = train_x, label = train_y)

      xgb_params_i <- list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        max_depth = as.integer(grid_i$max_depth),
        eta = as.numeric(grid_i$eta),
        gamma = as.numeric(grid_i$gamma),
        colsample_bytree = as.numeric(grid_i$colsample_bytree),
        min_child_weight = as.numeric(grid_i$min_child_weight),
        subsample = as.numeric(grid_i$subsample),
        nthread = 1
      )

      xgb_model_i <- xgboost::xgb.train(
        params = xgb_params_i,
        data = dtrain,
        nrounds = as.integer(grid_i$nrounds),
        verbose = 0
      )

      valid_prob <- predict(xgb_model_i, newdata = valid_x)

      xgb_fold_pred_list[[length(xgb_fold_pred_list) + 1]] <- tibble(
        fold_id = f,
        obs = valid_y,
        prob = as.numeric(valid_prob)
      )
    }

    xgb_cv_pred_i <- bind_rows(xgb_fold_pred_list)

    metric_i <- calculate_metrics(xgb_cv_pred_i$prob, xgb_cv_pred_i$obs)

    xgb_tuning_all[[length(xgb_tuning_all) + 1]] <- grid_i %>%
      mutate(
        pa_rep = r,
        TSS = metric_i$TSS[1],
        AUC = metric_i$AUC[1],
        Sensitivity = metric_i$sensitivity[1],
        Specificity = metric_i$specificity[1]
      )
  }
}

xgb_tuning_result <- bind_rows(xgb_tuning_all)

xgb_tuning_summary <- xgb_tuning_result %>%
  group_by(nrounds, max_depth, eta, gamma, colsample_bytree,
           min_child_weight, subsample) %>%
  summarise(
    TSS_mean = mean(TSS, na.rm = TRUE),
    TSS_sd = sd(TSS, na.rm = TRUE),
    AUC_mean = mean(AUC, na.rm = TRUE),
    Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
    Specificity_mean = mean(Specificity, na.rm = TRUE),
    n_tuning_reps = n_distinct(pa_rep),
    .groups = "drop"
  ) %>%
  arrange(desc(TSS_mean), desc(AUC_mean))

best_xgb_params <- xgb_tuning_summary %>%
  slice(1)

best_xgb_params

write.csv(
  xgb_tuning_result,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_XGBoost_tuning_raw.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  xgb_tuning_summary,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_XGBoost_tuning_summary.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  best_xgb_params,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_XGBoost_selected_parameters.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 12. LightGBM Hyperparameter Tuning ---------####

lgb_grid <- expand.grid(
  nrounds = c(50, 100),
  learning_rate = c(0.01, 0.05),
  num_leaves = c(3, 7),
  max_depth = c(2, 3),
  min_data_in_leaf = c(5, 10),
  feature_fraction = c(0.8),
  bagging_fraction = c(0.7),
  bagging_freq = c(1),
  lambda_l1 = c(0),
  lambda_l2 = c(0, 1),
  min_gain_to_split = c(0)
) %>%
  filter(num_leaves <= 2^max_depth)

lgb_tuning_all <- list()

for (r in seq_len(n_tuning_reps)) {
  message("LightGBM tuning with repeated pseudoabsence dataset: ", r)

  train_set_r <- ml_pa_dataset_list[[r]] %>%
    dplyr::select(presence, all_of(env_vars), fold_id)

  fold_values_lgb <- sort(unique(train_set_r$fold_id))

  for (g in seq_len(nrow(lgb_grid))) {
    message("  LightGBM grid ", g, " / ", nrow(lgb_grid))

    grid_i <- lgb_grid[g, ]
    lgb_fold_pred_list <- list()

    for (f in fold_values_lgb) {
      train_fold_i <- train_set_r %>% filter(fold_id != f)
      valid_fold_i <- train_set_r %>% filter(fold_id == f)

      train_y <- ifelse(train_fold_i$presence == "X1", 1, 0)
      valid_y <- ifelse(valid_fold_i$presence == "X1", 1, 0)

      train_x <- as.matrix(as.data.frame(train_fold_i[, env_vars]))
      valid_x <- as.matrix(as.data.frame(valid_fold_i[, env_vars]))

      dtrain <- lightgbm::lgb.Dataset(data = train_x, label = train_y)

      lgb_params_i <- list(
        objective = "binary",
        metric = "binary_logloss",
        learning_rate = as.numeric(grid_i$learning_rate),
        num_leaves = as.integer(grid_i$num_leaves),
        max_depth = as.integer(grid_i$max_depth),
        min_data_in_leaf = as.integer(grid_i$min_data_in_leaf),
        feature_fraction = as.numeric(grid_i$feature_fraction),
        bagging_fraction = as.numeric(grid_i$bagging_fraction),
        bagging_freq = as.integer(grid_i$bagging_freq),
        lambda_l1 = as.numeric(grid_i$lambda_l1),
        lambda_l2 = as.numeric(grid_i$lambda_l2),
        min_gain_to_split = as.numeric(grid_i$min_gain_to_split),
        num_threads = 1,
        verbosity = -1,
        force_col_wise = TRUE
      )

      lgb_model_i <- lightgbm::lgb.train(
        params = lgb_params_i,
        data = dtrain,
        nrounds = as.integer(grid_i$nrounds),
        verbose = -1
      )

      valid_prob <- predict(lgb_model_i, newdata = valid_x)

      lgb_fold_pred_list[[length(lgb_fold_pred_list) + 1]] <- tibble(
        fold_id = f,
        obs = valid_y,
        prob = as.numeric(valid_prob)
      )
    }

    lgb_cv_pred_i <- bind_rows(lgb_fold_pred_list)

    metric_i <- calculate_metrics(lgb_cv_pred_i$prob, lgb_cv_pred_i$obs)

    lgb_tuning_all[[length(lgb_tuning_all) + 1]] <- grid_i %>%
      mutate(
        pa_rep = r,
        TSS = metric_i$TSS[1],
        AUC = metric_i$AUC[1],
        Sensitivity = metric_i$sensitivity[1],
        Specificity = metric_i$specificity[1]
      )
  }
}

lgb_tuning_result <- bind_rows(lgb_tuning_all)

lgb_tuning_summary <- lgb_tuning_result %>%
  group_by(nrounds, learning_rate, num_leaves, max_depth, min_data_in_leaf,
                  feature_fraction, bagging_fraction, bagging_freq,
           lambda_l1, lambda_l2, min_gain_to_split) %>%
  summarise(
    TSS_mean = mean(TSS, na.rm = TRUE),
    TSS_sd = sd(TSS, na.rm = TRUE),
    AUC_mean = mean(AUC, na.rm = TRUE),
    Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
    Specificity_mean = mean(Specificity, na.rm = TRUE),
    n_tuning_reps = n_distinct(pa_rep),
    .groups = "drop"
  ) %>%
  arrange(desc(TSS_mean), desc(AUC_mean))

best_lgb_params <- lgb_tuning_summary %>%
  slice(1)

best_lgb_params

write.csv(
  lgb_tuning_result,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_LightGBM_tuning_raw.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  lgb_tuning_summary,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_LightGBM_tuning_summary.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  best_lgb_params,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_LightGBM_selected_parameters.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 13-1. RF Spatial Out-of-fold Predictions ---------####

# Fit with the selected hyperparameters.
fit_rf_selected_params <- function(train_df, env_vars, best_pars) {
  train_x <- train_df %>%
    dplyr::select(presence, all_of(env_vars)) %>%
    mutate(presence = factor(presence, levels = c("X0", "X1")))

  rf_model <- ranger::ranger(
    presence ~ .,
    data = train_x,
    probability = TRUE,
    num.trees = as.integer(best_pars$num_trees),
    mtry = as.integer(best_pars$mtry),
    splitrule = as.character(best_pars$splitrule),
    min.node.size = as.integer(best_pars$min.node.size),
    sample.fraction = as.numeric(best_pars$sample.fraction),
    importance = "permutation",
    seed = 777
  )

  return(rf_model)
}


# Training uses 3x pseudoabsences; validation uses the full held-out background.
rf_oof_list <- list()

set.seed(777)
for (r in seq_len(n_pa_reps)) {
  for (f in seq_len(n_folds)) {
    occ_train <- occ_folded %>% filter(fold_id != f)
    bg_train <- bg_folded %>% filter(fold_id != f)
    valid_set <- pa_model_data %>% filter(fold_id == f)

    n_bg <- min(nrow(bg_train), nrow(occ_train) * 3)
    bg_sample <- slice_sample(bg_train, n = n_bg, replace = FALSE)

    train_set <- bind_rows(
      occ_train %>% mutate(presence = factor("X1", levels = c("X0", "X1"))),
      bg_sample %>% mutate(presence = factor("X0", levels = c("X0", "X1")))
    )

    model_f <- fit_rf_selected_params(train_set, env_vars, best_rf_params)
    pred_input <- as.data.frame(valid_set[, env_vars])
    pred_f <- as.numeric(predict(model_f, data = pred_input)$predictions[, "X1"])

    rf_oof_list[[length(rf_oof_list) + 1]] <- valid_set %>%
      transmute(str_ID, fold_id, obs = as.integer(presence == "X1"),
                pa_rep = r, prob = pred_f)
  }
}

rf_oof_pred <- bind_rows(rf_oof_list) %>%
  group_by(str_ID, fold_id, obs) %>%
  summarise(prob_rf = mean(prob), .groups = "drop")

write.csv(
  rf_oof_pred,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_RF_spatial_OOF_predictions.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 13-2. GBM Spatial Out-of-fold Predictions ---------####

# Fit with the selected hyperparameters.
fit_gbm_selected_params <- function(train_df, env_vars, best_pars, bag_fraction = 0.70) {
  train_x <- train_df %>%
    mutate(y = ifelse(presence == "X1", 1, 0)) %>%
    dplyr::select(y, all_of(env_vars))

  gbm_model <- gbm::gbm(
    formula = y ~ .,
    data = train_x,
    distribution = "bernoulli",
    n.trees = as.integer(best_pars$n.trees),
    interaction.depth = as.integer(best_pars$interaction.depth),
    shrinkage = as.numeric(best_pars$shrinkage),
    n.minobsinnode = as.integer(best_pars$n.minobsinnode),
    bag.fraction = bag_fraction,
    train.fraction = 1.0,
    keep.data = FALSE,
    verbose = FALSE
  )

  attr(gbm_model, "n.trees") <- as.integer(best_pars$n.trees)
  return(gbm_model)
}


# Training uses 3x pseudoabsences; validation uses the full held-out background.
gbm_oof_list <- list()

set.seed(888)
for (r in seq_len(n_pa_reps)) {
  for (f in seq_len(n_folds)) {
    occ_train <- occ_folded %>% filter(fold_id != f)
    bg_train <- bg_folded %>% filter(fold_id != f)
    valid_set <- pa_model_data %>% filter(fold_id == f)

    n_bg <- min(nrow(bg_train), nrow(occ_train) * 3)
    bg_sample <- slice_sample(bg_train, n = n_bg, replace = FALSE)

    train_set <- bind_rows(
      occ_train %>% mutate(presence = factor("X1", levels = c("X0", "X1"))),
      bg_sample %>% mutate(presence = factor("X0", levels = c("X0", "X1")))
    )

    model_f <- fit_gbm_selected_params(train_set, env_vars, best_gbm_params)
    pred_input <- as.data.frame(valid_set[, env_vars])
    pred_f <- as.numeric(predict(model_f, newdata = pred_input,
                                   n.trees = attr(model_f, "n.trees"), type = "response"))

    gbm_oof_list[[length(gbm_oof_list) + 1]] <- valid_set %>%
      transmute(str_ID, fold_id, obs = as.integer(presence == "X1"),
                pa_rep = r, prob = pred_f)
  }
}

gbm_oof_pred <- bind_rows(gbm_oof_list) %>%
  group_by(str_ID, fold_id, obs) %>%
  summarise(prob_gbm = mean(prob), .groups = "drop")

write.csv(
  gbm_oof_pred,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_GBM_spatial_OOF_predictions.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 13-3. ANN Spatial Out-of-fold Predictions ---------####

# Fit with the selected hyperparameters.
fit_ann_selected_params <- function(train_df, env_vars, best_pars) {
  train_x <- train_df %>%
    dplyr::select(presence, all_of(env_vars)) %>%
    mutate(presence = factor(presence, levels = c("X0", "X1")))

  ann_grid_i <- best_pars %>%
    dplyr::select(size, decay) %>%
    mutate(
      size = as.integer(size),
      decay = as.numeric(decay)
    )

  ann_control <- caret::trainControl(
    method = "none",
    classProbs = TRUE
  )

  ann_model <- caret::train(
    presence ~ .,
    data = train_x,
    method = "nnet",
    trControl = ann_control,
    tuneGrid = ann_grid_i,
    preProcess = c("center", "scale"),
    trace = FALSE,
    maxit = 500,
    MaxNWts = 10000
  )

  return(ann_model)
}


# Training uses 3x pseudoabsences; validation uses the full held-out background.
ann_oof_list <- list()

set.seed(999)
for (r in seq_len(n_pa_reps)) {
  for (f in seq_len(n_folds)) {
    occ_train <- occ_folded %>% filter(fold_id != f)
    bg_train <- bg_folded %>% filter(fold_id != f)
    valid_set <- pa_model_data %>% filter(fold_id == f)

    n_bg <- min(nrow(bg_train), nrow(occ_train) * 3)
    bg_sample <- slice_sample(bg_train, n = n_bg, replace = FALSE)

    train_set <- bind_rows(
      occ_train %>% mutate(presence = factor("X1", levels = c("X0", "X1"))),
      bg_sample %>% mutate(presence = factor("X0", levels = c("X0", "X1")))
    )

    model_f <- fit_ann_selected_params(train_set, env_vars, best_ann_params)
    pred_input <- as.data.frame(valid_set[, env_vars])
    pred_f <- as.numeric(predict(model_f, newdata = pred_input, type = "prob")[, "X1"])

    ann_oof_list[[length(ann_oof_list) + 1]] <- valid_set %>%
      transmute(str_ID, fold_id, obs = as.integer(presence == "X1"),
                pa_rep = r, prob = pred_f)
  }
}

ann_oof_pred <- bind_rows(ann_oof_list) %>%
  group_by(str_ID, fold_id, obs) %>%
  summarise(prob_ann = mean(prob), .groups = "drop")

write.csv(
  ann_oof_pred,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ANN_spatial_OOF_predictions.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 13-4. XGBoost Spatial Out-of-fold Predictions ---------####

# Fit with the selected hyperparameters.
fit_xgb_selected_params <- function(train_df, env_vars, best_pars) {
  train_y <- ifelse(train_df$presence == "X1", 1, 0)
  train_x <- train_df %>%
    dplyr::select(all_of(env_vars)) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()

  dtrain <- xgboost::xgb.DMatrix(
    data = train_x,
    label = train_y
  )

  xgb_params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = as.integer(best_pars$max_depth[1]),
    eta = as.numeric(best_pars$eta[1]),
    gamma = as.numeric(best_pars$gamma[1]),
    colsample_bytree = as.numeric(best_pars$colsample_bytree[1]),
    min_child_weight = as.numeric(best_pars$min_child_weight[1]),
    subsample = as.numeric(best_pars$subsample[1]),
    nthread = 1
  )

  xgb_model <- xgboost::xgb.train(
    params = xgb_params,
    data = dtrain,
    nrounds = as.integer(best_pars$nrounds[1]),
    verbose = 0
  )

  attr(xgb_model, "best_pars") <- best_pars
  attr(xgb_model, "env_vars") <- env_vars
  return(xgb_model)
}


# Training uses 3x pseudoabsences; validation uses the full held-out background.
xgb_oof_list <- list()

set.seed(1001)
for (r in seq_len(n_pa_reps)) {
  for (f in seq_len(n_folds)) {
    occ_train <- occ_folded %>% filter(fold_id != f)
    bg_train <- bg_folded %>% filter(fold_id != f)
    valid_set <- pa_model_data %>% filter(fold_id == f)

    n_bg <- min(nrow(bg_train), nrow(occ_train) * 3)
    bg_sample <- slice_sample(bg_train, n = n_bg, replace = FALSE)

    train_set <- bind_rows(
      occ_train %>% mutate(presence = factor("X1", levels = c("X0", "X1"))),
      bg_sample %>% mutate(presence = factor("X0", levels = c("X0", "X1")))
    )

    model_f <- fit_xgb_selected_params(train_set, env_vars, best_xgb_params)
    pred_input <- as.data.frame(valid_set[, env_vars])
    pred_f <- as.numeric(predict(model_f, newdata = as.matrix(pred_input)))

    xgb_oof_list[[length(xgb_oof_list) + 1]] <- valid_set %>%
      transmute(str_ID, fold_id, obs = as.integer(presence == "X1"),
                pa_rep = r, prob = pred_f)
  }
}

xgb_oof_pred <- bind_rows(xgb_oof_list) %>%
  group_by(str_ID, fold_id, obs) %>%
  summarise(prob_xgb = mean(prob), .groups = "drop")

write.csv(
  xgb_oof_pred,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_XGBoost_spatial_OOF_predictions.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 13-5. LightGBM Spatial Out-of-fold Predictions ---------####

# Fit with the selected hyperparameters.
fit_lgb_selected_params <- function(train_df, env_vars, best_pars) {
  train_y <- ifelse(train_df$presence == "X1", 1, 0)
  train_x <- train_df %>%
    dplyr::select(all_of(env_vars)) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()

  dtrain <- lightgbm::lgb.Dataset(
    data = train_x,
    label = train_y
  )

  lgb_params <- list(
    objective = "binary",
    metric = "binary_logloss",
    learning_rate = as.numeric(best_pars$learning_rate[1]),
    num_leaves = as.integer(best_pars$num_leaves[1]),
    max_depth = as.integer(best_pars$max_depth[1]),
    min_data_in_leaf = as.integer(best_pars$min_data_in_leaf[1]),
    feature_fraction = as.numeric(best_pars$feature_fraction[1]),
    bagging_fraction = as.numeric(best_pars$bagging_fraction[1]),
    bagging_freq = as.integer(best_pars$bagging_freq[1]),
    lambda_l1 = as.numeric(best_pars$lambda_l1[1]),
    lambda_l2 = as.numeric(best_pars$lambda_l2[1]),
    min_gain_to_split = as.numeric(best_pars$min_gain_to_split[1]),
    seed = 777,
    num_threads = 1,
    verbosity = -1,
    force_col_wise = TRUE
  )

  lgb_model <- lightgbm::lgb.train(
    params = lgb_params,
    data = dtrain,
    nrounds = as.integer(best_pars$nrounds[1]),
    verbose = -1,
    serializable = TRUE
  )

  attr(lgb_model, "best_pars") <- best_pars
  attr(lgb_model, "env_vars") <- env_vars
  return(lgb_model)
}


# Training uses 3x pseudoabsences; validation uses the full held-out background.
lgb_oof_list <- list()

set.seed(1002)
for (r in seq_len(n_pa_reps)) {
  for (f in seq_len(n_folds)) {
    occ_train <- occ_folded %>% filter(fold_id != f)
    bg_train <- bg_folded %>% filter(fold_id != f)
    valid_set <- pa_model_data %>% filter(fold_id == f)

    n_bg <- min(nrow(bg_train), nrow(occ_train) * 3)
    bg_sample <- slice_sample(bg_train, n = n_bg, replace = FALSE)

    train_set <- bind_rows(
      occ_train %>% mutate(presence = factor("X1", levels = c("X0", "X1"))),
      bg_sample %>% mutate(presence = factor("X0", levels = c("X0", "X1")))
    )

    model_f <- fit_lgb_selected_params(train_set, env_vars, best_lgb_params)
    pred_input <- as.data.frame(valid_set[, env_vars])
    pred_f <- as.numeric(predict(model_f, newdata = as.matrix(pred_input), type = "response"))

    lgb_oof_list[[length(lgb_oof_list) + 1]] <- valid_set %>%
      transmute(str_ID, fold_id, obs = as.integer(presence == "X1"),
                pa_rep = r, prob = pred_f)
  }
}

lgb_oof_pred <- bind_rows(lgb_oof_list) %>%
  group_by(str_ID, fold_id, obs) %>%
  summarise(prob_lgb = mean(prob), .groups = "drop")

write.csv(
  lgb_oof_pred,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_LightGBM_spatial_OOF_predictions.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 14. Model Performance & Ensemble Weights ---------####

# Join predictions by stream segment.
ensemble_oof_predictions <- maxent_oof_pred %>%
  left_join(rf_oof_pred, by = c("str_ID", "fold_id", "obs")) %>%
  left_join(gbm_oof_pred, by = c("str_ID", "fold_id", "obs")) %>%
  left_join(ann_oof_pred, by = c("str_ID", "fold_id", "obs")) %>%
  left_join(xgb_oof_pred, by = c("str_ID", "fold_id", "obs")) %>%
  left_join(lgb_oof_pred, by = c("str_ID", "fold_id", "obs"))

prob_columns <- c("prob_maxent", "prob_rf", "prob_gbm", "prob_ann", "prob_xgb", "prob_lgb")
model_names <- c("MaxEnt", "RF", "GBM", "ANN", "XGBoost", "LightGBM")

# Pooled OOF metrics at the maxSSS threshold for each algorithm
model_metrics <- data.frame()

for (i in seq_along(prob_columns)) {
  metrics_i <- calculate_metrics(
    ensemble_oof_predictions[[prob_columns[i]]],
    ensemble_oof_predictions$obs
  )

  model_metrics <- bind_rows(
    model_metrics,
    mutate(metrics_i, model = model_names[i], .before = 1)
  )
}


# Normalize non-negative TSS values to sum to one.
weight_table <- model_metrics %>%
  dplyr::select(model, TSS) %>%
  mutate(weight_raw = ifelse(is.na(TSS) | TSS < 0, 0, TSS))

if (sum(weight_table$weight_raw) > 0) {
  weight_table$weight <- weight_table$weight_raw / sum(weight_table$weight_raw)
} else {
  weight_table$weight <- rep(1 / nrow(weight_table), nrow(weight_table))
}


# Ensemble OOF predictions and threshold
prob_matrix <- as.matrix(ensemble_oof_predictions[, prob_columns])
ensemble_oof_predictions$prob_ensemble <- as.numeric(prob_matrix %*% weight_table$weight)

ensemble_metrics <- calculate_metrics(
  ensemble_oof_predictions$prob_ensemble,
  ensemble_oof_predictions$obs
) %>%
  mutate(model = "Ensemble", .before = 1)

ensemble_maxsss_thr <- ensemble_metrics$maxsss_thr[1]

model_performance <- bind_rows(model_metrics, ensemble_metrics)

model_performance

write.csv(
  model_performance,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_model_performance.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  weight_table,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ensemble_weights.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  ensemble_oof_predictions,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ensemble_OOF_predictions.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 15-1. Final RF Models ---------####

# Fit ten models using independently sampled pseudoabsences.
rf_final_models <- vector("list", n_pa_reps)

set.seed(777)
for (r in seq_len(n_pa_reps)) {
  n_bg <- min(nrow(bg_folded), nrow(occ_folded) * 3)
  bg_sample <- slice_sample(bg_folded, n = n_bg, replace = FALSE)

  train_set <- bind_rows(
    occ_folded %>% mutate(presence = factor("X1", levels = c("X0", "X1"))),
    bg_sample %>% mutate(presence = factor("X0", levels = c("X0", "X1")))
  )

  rf_final_models[[r]] <- fit_rf_selected_params(
    train_set, env_vars, best_rf_params
  )
}

saveRDS(
  rf_final_models,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_RF_final_models.rds"
)



####--------- 15-2. Final GBM Models ---------####

# Fit ten models using independently sampled pseudoabsences.
gbm_final_models <- vector("list", n_pa_reps)

set.seed(888)
for (r in seq_len(n_pa_reps)) {
  n_bg <- min(nrow(bg_folded), nrow(occ_folded) * 3)
  bg_sample <- slice_sample(bg_folded, n = n_bg, replace = FALSE)

  train_set <- bind_rows(
    occ_folded %>% mutate(presence = factor("X1", levels = c("X0", "X1"))),
    bg_sample %>% mutate(presence = factor("X0", levels = c("X0", "X1")))
  )

  gbm_final_models[[r]] <- fit_gbm_selected_params(
    train_set, env_vars, best_gbm_params
  )
}

saveRDS(
  gbm_final_models,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_GBM_final_models.rds"
)



####--------- 15-3. Final ANN Models ---------####

# Fit ten models using independently sampled pseudoabsences.
ann_final_models <- vector("list", n_pa_reps)

set.seed(999)
for (r in seq_len(n_pa_reps)) {
  n_bg <- min(nrow(bg_folded), nrow(occ_folded) * 3)
  bg_sample <- slice_sample(bg_folded, n = n_bg, replace = FALSE)

  train_set <- bind_rows(
    occ_folded %>% mutate(presence = factor("X1", levels = c("X0", "X1"))),
    bg_sample %>% mutate(presence = factor("X0", levels = c("X0", "X1")))
  )

  ann_final_models[[r]] <- fit_ann_selected_params(
    train_set, env_vars, best_ann_params
  )
}

saveRDS(
  ann_final_models,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ANN_final_models.rds"
)



####--------- 15-4. Final XGBoost Models ---------####

# Fit ten models using independently sampled pseudoabsences.
xgb_final_models <- vector("list", n_pa_reps)

set.seed(1001)
for (r in seq_len(n_pa_reps)) {
  n_bg <- min(nrow(bg_folded), nrow(occ_folded) * 3)
  bg_sample <- slice_sample(bg_folded, n = n_bg, replace = FALSE)

  train_set <- bind_rows(
    occ_folded %>% mutate(presence = factor("X1", levels = c("X0", "X1"))),
    bg_sample %>% mutate(presence = factor("X0", levels = c("X0", "X1")))
  )

  xgb_final_models[[r]] <- fit_xgb_selected_params(
    train_set, env_vars, best_xgb_params
  )
}

for (r in seq_along(xgb_final_models)) {
  xgboost::xgb.save(
    xgb_final_models[[r]],
    fname = paste0("./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_XGBoost_rep_", r, ".json")
  )
}



####--------- 15-5. Final LightGBM Models ---------####

# Fit ten models using independently sampled pseudoabsences.
lgb_final_models <- vector("list", n_pa_reps)

set.seed(1002)
for (r in seq_len(n_pa_reps)) {
  n_bg <- min(nrow(bg_folded), nrow(occ_folded) * 3)
  bg_sample <- slice_sample(bg_folded, n = n_bg, replace = FALSE)

  train_set <- bind_rows(
    occ_folded %>% mutate(presence = factor("X1", levels = c("X0", "X1"))),
    bg_sample %>% mutate(presence = factor("X0", levels = c("X0", "X1")))
  )

  lgb_final_models[[r]] <- fit_lgb_selected_params(
    train_set, env_vars, best_lgb_params
  )
}

for (r in seq_along(lgb_final_models)) {
  lightgbm::lgb.save(
    lgb_final_models[[r]],
    filename = paste0("./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_LightGBM_rep_", r, ".txt")
  )
}



####--------- 16. Current Distribution Prediction ---------####

# Average repeated predictions within each algorithm, then apply ensemble weights.
predict_ensemble <- function(pred_input) {
  pred_input <- as.data.frame(pred_input[, env_vars, drop = FALSE])

  pred_maxent <- predict(mod_maxent, x = pred_input, args = "outputformat=cloglog")

  pred_rf <- lapply(rf_final_models, function(model) {
    as.numeric(predict(model, data = pred_input)$predictions[, "X1"])
  })

  pred_gbm <- lapply(gbm_final_models, function(model) {
    as.numeric(predict(model, newdata = pred_input,
                       n.trees = attr(model, "n.trees"), type = "response"))
  })

  pred_ann <- lapply(ann_final_models, function(model) {
    as.numeric(predict(model, newdata = pred_input, type = "prob")[, "X1"])
  })

  pred_xgb <- lapply(xgb_final_models, function(model) {
    as.numeric(predict(model, newdata = as.matrix(pred_input)))
  })

  pred_lgb <- lapply(lgb_final_models, function(model) {
    as.numeric(predict(model, newdata = as.matrix(pred_input), type = "response"))
  })

  predictions <- data.frame(
    prob_maxent = as.numeric(pred_maxent),
    prob_rf = rowMeans(do.call(cbind, pred_rf), na.rm = TRUE),
    prob_gbm = rowMeans(do.call(cbind, pred_gbm), na.rm = TRUE),
    prob_ann = rowMeans(do.call(cbind, pred_ann), na.rm = TRUE),
    prob_xgb = rowMeans(do.call(cbind, pred_xgb), na.rm = TRUE),
    prob_lgb = rowMeans(do.call(cbind, pred_lgb), na.rm = TRUE)
  )

  predictions$prob_ensemble <- as.numeric(as.matrix(predictions) %*% weight_table$weight)

  return(predictions)
}


# Stream lengths use the supplied length attribute in metres.
network_index <- streamnetwork %>%
  st_drop_geometry() %>%
  transmute(str_ID, BBSNCD, SBSNCD,
            stream_length_m = as.numeric(.data$length),
            in_accessible_area = SBSNCD %in% accessible_basins)

proj_with_id <- streamnetwork2 %>%
  st_drop_geometry() %>%
  dplyr::select(str_ID, longitude = mean_lon, latitude = mean_lat, all_of(env_vars)) %>%
  mutate(across(all_of(env_vars), as.numeric)) %>%
  na.omit() %>%
  filter(if_all(all_of(env_vars), is.finite)) %>%
  distinct(str_ID, .keep_all = TRUE) %>%
  as.data.frame()

pred_input <- proj_with_id[, env_vars, drop = FALSE]
pred_present <- predict_ensemble(pred_input)

pred_present <- bind_cols(
  proj_with_id[, c("str_ID", "longitude", "latitude")],
  pred_present
) %>%
  mutate(suitable = as.integer(prob_ensemble >= ensemble_maxsss_thr))

# NA indicates an unmodeled segment, not unsuitable habitat.
prediction_present <- network_index %>%
  left_join(pred_present, by = "str_ID") %>%
  mutate(scenario_id = "present", ssp_scenario = "present", time = "2000-2019",
         threshold = ensemble_maxsss_thr, .before = 1)

write.csv(
  prediction_present,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_current_predictions.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 17. Future Distribution Predictions ---------####

# Non-climatic predictors remain unchanged.
base_vars <- c("slp_mean", "slp_CV", "ws_area", "forest", "weir_dens")

# BIO05 and BIO13 columns for each SSP and period
scenarios <- list(
  s1_21_40 = c("s1_24_b05", "s1_24_b13"),
  s1_41_60 = c("s1_46_b05", "s1_46_b13"),
  s1_61_80 = c("s1_68_b05", "s1_68_b13"),
  s1_81_100 = c("s1_81_b05", "s1_81_b13"),

  s2_21_40 = c("s2_24_b05", "s2_24_b13"),
  s2_41_60 = c("s2_46_b05", "s2_46_b13"),
  s2_61_80 = c("s2_68_b05", "s2_68_b13"),
  s2_81_100 = c("s2_81_b05", "s2_81_b13"),

  s5_21_40 = c("s5_24_b05", "s5_24_b13"),
  s5_41_60 = c("s5_46_b05", "s5_46_b13"),
  s5_61_80 = c("s5_68_b05", "s5_68_b13"),
  s5_81_100 = c("s5_81_b05", "s5_81_b13")
)

ssp_levels <- c("SSP126", "SSP245", "SSP585")
ssp_labels <- c(SSP126 = "SSP1-2.6", SSP245 = "SSP2-4.5", SSP585 = "SSP5-8.5")
period_levels <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")

scenario_info <- data.frame(
  scenario_id = names(scenarios),
  ssp_scenario = rep(ssp_levels, each = 4),
  time = rep(period_levels, 3)
)

prediction_results <- list(present = prediction_present)

for (i in seq_along(scenarios)) {
  scenario_name <- names(scenarios)[i]
  climate_columns <- scenarios[[i]]

  proj_future <- streamnetwork2 %>%
    st_drop_geometry() %>%
    dplyr::select(str_ID, longitude = mean_lon, latitude = mean_lat,
                  all_of(base_vars),
                  p_bio05 = all_of(climate_columns[1]),
                  p_bio13 = all_of(climate_columns[2])) %>%
    mutate(across(all_of(env_vars), as.numeric)) %>%
    na.omit() %>%
    filter(if_all(all_of(env_vars), is.finite)) %>%
    distinct(str_ID, .keep_all = TRUE) %>%
    as.data.frame()

  pred_future <- predict_ensemble(proj_future[, env_vars, drop = FALSE])

  pred_future <- bind_cols(
    proj_future[, c("str_ID", "longitude", "latitude")],
    pred_future
  ) %>%
    mutate(suitable = as.integer(prob_ensemble >= ensemble_maxsss_thr))

  prediction_results[[scenario_name]] <- network_index %>%
    left_join(pred_future, by = "str_ID") %>%
    mutate(scenario_id = scenario_name,
           ssp_scenario = scenario_info$ssp_scenario[i],
           time = scenario_info$time[i],
           threshold = ensemble_maxsss_thr,
           .before = 1)
}

prediction_long <- bind_rows(prediction_results)

write.csv(
  prediction_long,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ensemble_predictions_long.csv",
  row.names = FALSE,
  na = "NA"
)


####--------- 17-1. Export Stream-network Shapefile ---------####

# Keep identifiers, basin information, stream order, length, and coordinates.
# Retain all rows in streamnetwork, including segments outside the accessible basins.
aux_streamnetwork <- streamnetwork %>%
  dplyr::select(
    str_ID,
    any_of(c("BBSNCD", "BBSNNM", "MBSNCD", "MBSNNM", "SBSNCD", "SBSNNM",
             "st_order", "length", "mean_lon", "mean_lat"))
  )


# OM_p: observed presence, not present-day predicted suitability.
# 1 = occurrence record; 0 = no occurrence record in a currently modeled segment.
# no_model = no current prediction and no occurrence record.
present_shp_value <- prediction_present$suitable[
  match(aux_streamnetwork$str_ID, prediction_present$str_ID)
]

aux_streamnetwork <- aux_streamnetwork %>%
  mutate(
    OM_p = case_when(
      str_ID %in% occ$str_ID ~ "1",
      is.na(present_shp_value) ~ "no_model",
      TRUE ~ "0"
    )
  )

# To use present-day predicted suitability for OM_p instead, replace it with:
# aux_streamnetwork$OM_p <- coalesce(as.character(present_shp_value), "no_model")


# Future fields: s1/s2/s5 = SSP1-2.6/SSP2-4.5/SSP5-8.5.
# 24/46/68/81 = 2021-2040/2041-2060/2061-2080/2081-2100.
scenario_columns <- c(
  OM_s1_24 = "s1_21_40",
  OM_s1_46 = "s1_41_60",
  OM_s1_68 = "s1_61_80",
  OM_s1_81 = "s1_81_100",
  OM_s2_24 = "s2_21_40",
  OM_s2_46 = "s2_41_60",
  OM_s2_68 = "s2_61_80",
  OM_s2_81 = "s2_81_100",
  OM_s5_24 = "s5_21_40",
  OM_s5_46 = "s5_41_60",
  OM_s5_68 = "s5_61_80",
  OM_s5_81 = "s5_81_100"
)

for (column_name in names(scenario_columns)) {
  scenario_result <- prediction_results[[scenario_columns[[column_name]]]]
  scenario_value <- scenario_result$suitable[
    match(aux_streamnetwork$str_ID, scenario_result$str_ID)
  ]

  # Store all three categories as text in the shapefile.
  aux_streamnetwork[[column_name]] <- coalesce(as.character(scenario_value), "no_model")
}


# Export basic attributes plus 13 species fields; retain the original geometry.
sf::st_write(
  aux_streamnetwork,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_streamnetwork_status.shp",
  driver = "ESRI Shapefile",
  layer_options = "ENCODING=UTF-8",
  delete_layer = TRUE,
  quiet = TRUE
)

write.csv(
  aux_streamnetwork %>% st_drop_geometry(),
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_streamnetwork_status.csv",
  row.names = FALSE,
  na = "NA",
  fileEncoding = "UTF-8"
)



####--------- 18. Current & Future Binary Maps ---------####

# Stream width by stream order
stream_order_col <- intersect(
  c("streamorde", "stream_order", "streamorde_", "str_order", "order", "ord_strah", "strahler"),
  names(map_streamnetwork)
)[1]

if (is.na(stream_order_col)) {
  map_streamnetwork <- map_streamnetwork %>% mutate(plot_stream_order = 1)
} else {
  map_streamnetwork <- map_streamnetwork %>%
    mutate(plot_stream_order = as.numeric(.data[[stream_order_col]]))
}

map_streamnetwork <- map_streamnetwork %>%
  mutate(
    plot_stream_order = if_else(is.na(plot_stream_order), 1, plot_stream_order),
    stream_linewidth = case_when(
      plot_stream_order <= 1 ~ 0.12,
      plot_stream_order == 2 ~ 0.18,
      plot_stream_order == 3 ~ 0.26,
      plot_stream_order == 4 ~ 0.36,
      plot_stream_order == 5 ~ 0.48,
      plot_stream_order == 6 ~ 0.62,
      TRUE ~ 0.76
    )
  )

# Crop extent based on stream segments used in the analysis area.
crop_stream <- map_streamnetwork %>%
  filter(str_ID %in% streamnetwork2$str_ID)

crop_bbox <- st_bbox(crop_stream)
crop_xpad <- as.numeric(crop_bbox["xmax"] - crop_bbox["xmin"]) * 0.03
crop_ypad <- as.numeric(crop_bbox["ymax"] - crop_bbox["ymin"]) * 0.03
crop_xlim <- c(crop_bbox["xmin"] - crop_xpad, crop_bbox["xmax"] + crop_xpad)
crop_ylim <- c(crop_bbox["ymin"] - crop_ypad, crop_bbox["ymax"] + crop_ypad)

# Present map layers
present_status <- prediction_present %>%
  st_drop_geometry() %>%
  dplyr::select(str_ID, suitable)

present_map <- map_streamnetwork %>%
  left_join(present_status, by = "str_ID")

present_unsuitable <- present_map %>% filter(suitable == 0)
present_suitable <- present_map %>% filter(suitable == 1)

present_binary_fig <- ggplot() +
  geom_sf(data = crop_stream, colour = "gray85", aes(linewidth = stream_linewidth), show.legend = FALSE) +
  geom_sf(data = present_unsuitable, colour = "#88deeb", aes(linewidth = stream_linewidth), show.legend = FALSE) +
  geom_sf(data = present_suitable, colour = "#d7191c", aes(linewidth = stream_linewidth), show.legend = FALSE) +
  geom_sf(data = map_watershed, fill = NA, colour = "gray30", linewidth = 0.45) +
  scale_linewidth_identity() +
  coord_sf(xlim = crop_xlim, ylim = crop_ylim, expand = FALSE) +
  ggspatial::annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.65,
    line_width = 0.35,
    pad_x = grid::unit(0.18, "cm"),
    pad_y = grid::unit(0.18, "cm")
  ) +
  ggspatial::annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = grid::unit(0.9, "cm"),
    width = grid::unit(0.9, "cm"),
    pad_x = grid::unit(0.18, "cm"),
    pad_y = grid::unit(0.18, "cm"),
    style = ggspatial::north_arrow_orienteering
  ) +
  labs(title = "Present") +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.background = element_rect(fill = "white", colour = "gray30", linewidth = 0.8),
    plot.margin = margin(5, 5, 5, 5)
  )

# Future maps in a cropped 3 x 4 panel
future_prediction_long <- bind_rows(prediction_results[-1]) %>%
  mutate(
    ssp_scenario = factor(ssp_scenario,
                          levels = c("SSP126", "SSP245", "SSP585"),
                          labels = c("SSP1-2.6", "SSP2-4.5", "SSP5-8.5")),
    time = factor(time,
                  levels = c("2021-2040", "2041-2060", "2061-2080", "2081-2100"))
  )

future_status <- future_prediction_long %>%
  st_drop_geometry() %>%
  dplyr::select(str_ID, ssp_scenario, time, suitable)

future_map <- map_streamnetwork %>%
  left_join(future_status, by = "str_ID") %>%
  mutate(prediction_class = case_when(
    suitable == 1 ~ "Suitable",
    suitable == 0 ~ "Unsuitable",
    TRUE ~ NA_character_
  ))

future_binary_fig <- ggplot() +
  geom_sf(data = crop_stream, colour = "gray85", aes(linewidth = stream_linewidth), show.legend = FALSE) +
  geom_sf(data = future_map %>% filter(prediction_class == "Unsuitable"),
          colour = "#88deeb", aes(linewidth = stream_linewidth), show.legend = FALSE) +
  geom_sf(data = future_map %>% filter(prediction_class == "Suitable"),
          colour = "#d7191c", aes(linewidth = stream_linewidth), show.legend = FALSE) +
  geom_sf(data = map_watershed, fill = NA, colour = "gray30", linewidth = 0.40) +
  scale_linewidth_identity() +
  coord_sf(xlim = crop_xlim, ylim = crop_ylim, expand = FALSE) +
  facet_grid(ssp_scenario ~ time, switch = "y") +
  ggspatial::annotation_scale(
    location = "bl",
    width_hint = 0.24,
    text_cex = 0.42,
    line_width = 0.24,
    pad_x = grid::unit(0.08, "cm"),
    pad_y = grid::unit(0.08, "cm")
  ) +
  ggspatial::annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = grid::unit(0.48, "cm"),
    width = grid::unit(0.48, "cm"),
    pad_x = grid::unit(0.08, "cm"),
    pad_y = grid::unit(0.08, "cm"),
    style = ggspatial::north_arrow_orienteering
  ) +
  theme_void(base_size = 9) +
  theme(
    panel.spacing = grid::unit(0, "lines"),
    strip.background = element_rect(fill = "gray80", colour = "gray30", linewidth = 0.8),
    strip.placement = "outside",
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y.right = element_text(face = "bold", size = 10, angle = -90),
    panel.border = element_rect(fill = NA, colour = "gray30", linewidth = 0.45),
    plot.background = element_rect(fill = "white", colour = NA)
  )

legend_binary_fig <- ggplot(
  data.frame(class = factor(c("Suitable", "Unsuitable"), levels = c("Suitable", "Unsuitable")), x = 1:2, y = 1)
) +
  geom_tile(aes(x = x, y = y, fill = class), width = 0.6, height = 0.4, colour = "gray30") +
  scale_fill_manual(values = c(Suitable = "#d7191c", Unsuitable = "#88deeb"), guide = "none") +
  geom_text(aes(x = x + 0.35, y = y, label = class), hjust = 0, fontface = "bold", size = 5) +
  coord_cartesian(xlim = c(0.6, 3.3), ylim = c(0.7, 1.3), clip = "off") +
  theme_void()

upper_binary_panel <- ggarrange(
  present_binary_fig,
  future_binary_fig,
  ncol = 2,
  widths = c(1.25, 1.90),
  align = "h"
)

binary_map_panel <- ggarrange(
  upper_binary_panel,
  legend_binary_fig,
  ncol = 1,
  heights = c(1, 0.08)
)

binary_map_panel

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_current_future_binary_maps.png",
  plot = binary_map_panel,
  width = 12.5,
  height = 7.4,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_current_future_binary_maps.svg",
  plot = binary_map_panel,
  width = 12.5,
  height = 7.4,
  device = grDevices::svg,
  bg = "white"
)


####--------- 19. Suitable Stream Length & Habitat Change ---------####

# Summarize only stream segments with valid predictions.
habitat_summary <- prediction_long %>%
  group_by(scenario_id, ssp_scenario, time) %>%
  summarise(
    n_total_projected = sum(!is.na(suitable)),
    n_suitable = sum(suitable == 1, na.rm = TRUE),
    n_unsuitable = sum(suitable == 0, na.rm = TRUE),
    n_missing_in_accessible_area = sum(in_accessible_area & is.na(suitable)),
    projected_total_length_km = sum(stream_length_m[!is.na(suitable)]) / 1000,
    suitable_length_km = sum(stream_length_m[which(suitable == 1)]) / 1000,
    .groups = "drop"
  ) %>%
  mutate(suitable_length_pct_of_projected =
           100 * suitable_length_km / projected_total_length_km)

present_summary <- habitat_summary %>% filter(scenario_id == "present")
current_suitable_km <- present_summary$suitable_length_km[1]

# Positive percentage_decrease indicates habitat loss.
habitat_summary <- habitat_summary %>%
  mutate(
    current_suitable_length_km = current_suitable_km,
    change_suitable_length_km = suitable_length_km - current_suitable_km,
    percentage_decrease = if (current_suitable_km > 0)
      100 * (current_suitable_km - suitable_length_km) / current_suitable_km else NA_real_
  ) %>%
  arrange(match(scenario_id, c("present", scenario_info$scenario_id)))

habitat_summary

write.csv(
  habitat_summary,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_habitat_change_summary.csv",
  row.names = FALSE,
  na = "NA"
)


# Percentage of valid projected length that is suitable
plot_current <- present_summary[rep(1, length(ssp_levels)), ]
plot_current$ssp_scenario <- ssp_levels
plot_current$period <- "Present"

plot_future <- habitat_summary %>%
  filter(scenario_id != "present") %>%
  mutate(period = time)

length_plot_data <- bind_rows(plot_current, plot_future) %>%
  mutate(period = factor(period, levels = c("Present", period_levels)),
         ssp_scenario = factor(ssp_scenario, levels = ssp_levels))

length_plot <- ggplot(
  length_plot_data,
  aes(period, suitable_length_pct_of_projected, colour = ssp_scenario, group = ssp_scenario)
) +
  geom_line(linewidth = 0.9) +
  geom_point(data = filter(length_plot_data, period != "Present"), size = 2.6) +
  geom_point(data = length_plot_data[1, ], colour = "black",
             shape = 21, fill = "white", size = 3) +
  scale_colour_manual(
    values = c(SSP126 = "#00A6A6", SSP245 = "#F28E1C", SSP585 = "#A50000"),
    labels = ssp_labels,
    name = "Scenario"
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = NULL, y = "Suitable stream length (% of predicted length)",
       title = "Oncorhynchus masou masou") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "italic", hjust = 0.5)
  )

length_plot

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_suitable_stream_length.png",
  plot = length_plot,
  width = 7.2,
  height = 5.4,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_suitable_stream_length.svg",
  plot = length_plot,
  width = 7.2,
  height = 5.4,
  device = grDevices::svg,
  bg = "white"
)



####--------- 20. Ensemble Permutation Importance ---------####

# Evaluate final models on current presence-background data, not OOF predictions.
# Use the same apparent maxSSS threshold for all permutations.
n_permutations <- 30

baseline_x <- as.data.frame(pa_model_data[, env_vars])
baseline_obs <- as.integer(pa_model_data$presence == "X1")
baseline_prob <- predict_ensemble(baseline_x)$prob_ensemble
baseline_metrics <- calculate_metrics(baseline_prob, baseline_obs)

baseline_auc <- baseline_metrics$AUC[1]
baseline_tss <- baseline_metrics$TSS[1]
baseline_thr <- baseline_metrics$maxsss_thr[1]


# One row per predictor and permutation
importance_raw <- expand.grid(
  permutation_rep = seq_len(n_permutations),
  variable = env_vars,
  stringsAsFactors = FALSE
)
importance_raw$AUC_drop <- NA_real_
importance_raw$TSS_drop <- NA_real_

set.seed(888)
for (i in seq_len(nrow(importance_raw))) {
  variable_i <- importance_raw$variable[i]
  permuted_x <- baseline_x
  permuted_x[[variable_i]] <- sample(permuted_x[[variable_i]], replace = FALSE)

  permuted_prob <- predict_ensemble(permuted_x)$prob_ensemble
  permuted_metrics <- calculate_metrics(permuted_prob, baseline_obs, baseline_thr)

  importance_raw$AUC_drop[i] <- baseline_auc - permuted_metrics$AUC[1]
  importance_raw$TSS_drop[i] <- baseline_tss - permuted_metrics$TSS[1]
}

ensemble_importance_summary <- importance_raw %>%
  group_by(variable) %>%
  summarise(
    AUC_drop_mean = mean(AUC_drop),
    AUC_drop_sd = sd(AUC_drop),
    TSS_drop_mean = mean(TSS_drop),
    TSS_drop_sd = sd(TSS_drop),
    n_permutations = n(),
    .groups = "drop"
  )

auc_total <- sum(pmax(ensemble_importance_summary$AUC_drop_mean, 0))
tss_total <- sum(pmax(ensemble_importance_summary$TSS_drop_mean, 0))

ensemble_importance_summary <- ensemble_importance_summary %>%
  mutate(
    AUC_importance_pct = if (auc_total > 0)
      100 * pmax(AUC_drop_mean, 0) / auc_total else NA_real_,
    TSS_importance_pct = if (tss_total > 0)
      100 * pmax(TSS_drop_mean, 0) / tss_total else NA_real_
  ) %>%
  arrange(desc(AUC_drop_mean))

ensemble_importance_summary

write.csv(
  importance_raw,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ensemble_permutation_importance_raw.csv",
  row.names = FALSE,
  na = "NA"
)

write.csv(
  ensemble_importance_summary,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ensemble_permutation_importance_summary.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 21-1. Accumulated Local Effects ---------####

ale_n_bins <- 15
ale_n_boot <- 100
ale_top_n <- 6


# Quantile intervals; ties may reduce the number of bins.
make_ale_breaks <- function(x, n_bins = 20) {
  x <- as.numeric(x)
  x_finite <- x[is.finite(x)]

  if (length(unique(x_finite)) < 3) {
    return(NULL)
  }

  breaks <- stats::quantile(
    x_finite,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE,
    names = FALSE,
    type = 8
  )

  breaks <- sort(unique(as.numeric(breaks)))

  if (length(breaks) < 3) {
    breaks <- pretty(x_finite, n = min(n_bins, 10))
    breaks <- sort(unique(as.numeric(breaks)))
  }

  if (length(breaks) < 3) {
    return(NULL)
  }

  breaks[1] <- min(x_finite, na.rm = TRUE)
  breaks[length(breaks)] <- max(x_finite, na.rm = TRUE)
  breaks <- sort(unique(breaks))

  if (length(breaks) < 3) {
    return(NULL)
  }

  return(breaks)
}


# Calculate and center one-dimensional ALE for the fitted ensemble.
calculate_ale_1d <- function(data, variable, breaks) {
  data <- data %>%
    as.data.frame() %>%
    dplyr::select(all_of(env_vars)) %>%
    mutate(across(everything(), as.numeric))

  bin_id <- cut(data[[variable]], breaks = breaks,
                include.lowest = TRUE, labels = FALSE)

  ale_result <- data.frame(
    variable = variable,
    bin = seq_len(length(breaks) - 1),
    x_lower = head(breaks, -1),
    x_upper = tail(breaks, -1),
    x_mid = (head(breaks, -1) + tail(breaks, -1)) / 2,
    n = 0,
    local_effect = NA_real_
  )

  for (k in seq_len(nrow(ale_result))) {
    idx <- which(bin_id == k)
    if (length(idx) == 0) next

    data_lower <- data[idx, env_vars, drop = FALSE]
    data_upper <- data[idx, env_vars, drop = FALSE]
    data_lower[[variable]] <- breaks[k]
    data_upper[[variable]] <- breaks[k + 1]

    pred_lower <- predict_ensemble(data_lower)$prob_ensemble
    pred_upper <- predict_ensemble(data_upper)$prob_ensemble

    ale_result$n[k] <- length(idx)
    ale_result$local_effect[k] <- mean(pred_upper - pred_lower, na.rm = TRUE)
  }

  ale_result <- ale_result %>%
    filter(n > 0, !is.na(local_effect)) %>%
    arrange(bin) %>%
    mutate(ale_raw = cumsum(local_effect))

  ale_center <- weighted.mean(ale_result$ale_raw, w = ale_result$n, na.rm = TRUE)

  ale_result %>% mutate(ale_centered = ale_raw - ale_center)
}


# Select six predictors using AUC-based permutation importance.
ale_input_data <- as.data.frame(proj_with_id[, env_vars, drop = FALSE])

ale_top_vars <- ensemble_importance_summary %>%
  filter(variable %in% env_vars, is.finite(AUC_drop_mean)) %>%
  arrange(desc(AUC_drop_mean), desc(TSS_drop_mean)) %>%
  slice_head(n = ale_top_n) %>%
  pull(variable)

ale_results <- list()
ale_breaks <- list()

for (v in ale_top_vars) {
  breaks_v <- make_ale_breaks(ale_input_data[[v]], n_bins = ale_n_bins)
  if (is.null(breaks_v)) next

  ale_results[[v]] <- calculate_ale_1d(ale_input_data, v, breaks_v)
  ale_breaks[[v]] <- breaks_v
}

ale_result_long <- bind_rows(ale_results)

write.csv(
  ale_result_long,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ensemble_ALE_top6_variables.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 21-2. ALE Bootstrap Intervals ---------####

# Resample stream segments; keep fitted models and breakpoints unchanged.
set.seed(20260602)
ale_boot_list <- list()

for (v in names(ale_breaks)) {
  for (b in seq_len(ale_n_boot)) {
    idx <- sample(seq_len(nrow(ale_input_data)),
                  size = nrow(ale_input_data), replace = TRUE)

    ale_boot_list[[length(ale_boot_list) + 1]] <- calculate_ale_1d(
      ale_input_data[idx, , drop = FALSE],
      v,
      ale_breaks[[v]]
    ) %>%
      mutate(bootstrap_rep = b)
  }
}

ale_boot_long <- bind_rows(ale_boot_list)

ale_ci <- ale_boot_long %>%
  group_by(variable, bin, x_mid) %>%
  summarise(
    ale_centered_lwr = quantile(ale_centered, 0.025, na.rm = TRUE),
    ale_centered_upr = quantile(ale_centered, 0.975, na.rm = TRUE),
    n_boot_valid = sum(is.finite(ale_centered)),
    .groups = "drop"
  )

ale_result_ci <- ale_result_long %>%
  left_join(ale_ci, by = c("variable", "bin", "x_mid"))

write.csv(
  ale_result_ci,
  "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ensemble_ALE_top6_variables_with_bootstrap_CI.csv",
  row.names = FALSE,
  na = "NA"
)



####--------- 21-3. ALE Figure ---------####

# Show forest cover as a percentage; retain the model scale in CSV outputs.
forest_scale <- if (all(ale_input_data$forest >= 0 & ale_input_data$forest <= 1)) 100 else 1

predictor_labels <- c(
  slp_mean = "Mean stream slope (degrees)",
  slp_CV = "Slope CV",
  ws_area = "Upstream watershed area (km2)",
  forest = "Forest cover (%)",
  weir_dens = "Weir density (per km)",
  p_bio05 = "BIO05 (degrees C)",
  p_bio13 = "BIO13 (mm)"
)

ale_plot_data <- ale_result_ci %>%
  mutate(x_display = ifelse(variable == "forest", x_mid * forest_scale, x_mid),
         variable = factor(variable, levels = ale_top_vars))

ale_plot <- ggplot(ale_plot_data, aes(x_display, ale_centered)) +
  geom_ribbon(aes(ymin = ale_centered_lwr, ymax = ale_centered_upr), alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ variable, ncol = 3, scales = "free_x",
              labeller = as_labeller(predictor_labels)) +
  labs(x = NULL, y = "ALE", title = "Oncorhynchus masou masou") +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "italic", hjust = 0.5),
    legend.position = "none"
  )

ale_plot

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ensemble_ALE_top6.png",
  plot = ale_plot,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "./data_outcome/O_masou_6model_ensemble_V2.10/O_masou_ensemble_ALE_top6.svg",
  plot = ale_plot,
  width = 11,
  height = 7,
  device = grDevices::svg,
  bg = "white"
)
