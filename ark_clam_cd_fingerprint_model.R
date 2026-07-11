# ============================================================
# 0. Packages and input data
# ============================================================
library(readxl)
library(tidyverse)
library(glmnet)
library(pracma)
library(waveslim)
library(cowplot)

has_pracma <- TRUE
has_waveslim <- TRUE

set.seed(123)

protein_raw <- read_excel(
  "ark_clam_hplc_uv_profiles.xlsx",
  sheet = "Sheet1"
)

icpms_raw <- read_excel(
  "ark_clam_sec_icp_ms_profiles.xlsx",
  sheet = "Sheet1"
)

metal_raw <- read_excel(
  "ark_clam_subcellular_metal_concentrations.xlsx",
  sheet = "Sheet1"
)
# HPLC/SEC settings
bin_width <- 1
max_time  <- 60
n_quantile_bins <- 30

# Main target
target_metal <- "Cd"
target_col   <- "Cd_total_ug_g"
target_icpms_signal <- "Cd113"

# Metal concentration columns in "cytosol and particle.xlsx"
metal_total_map <- tribble(
  ~metal, ~back_col,         ~accumulated_col,         ~total_col,
  "Ni",   "Ni60_back_ug_g",  "accumulated_Ni61_ug_g",  "Ni_total_ug_g",
  "Cu",   "Cu63_back_ug_g",  "accumulated_Cu65_ug_g",  "Cu_total_ug_g",
  "Zn",   "Zn66_back_ug_g",  "accumulated_Zn68_ug_g",  "Zn_total_ug_g",
  "Cd",   "Cd111_back_ug_g", "accumulated_Cd113_ug_g", "Cd_total_ug_g",
  "Pb",   "Pb207_back_ug_g", "accumulated_Pb206_ug_g", "Pb_total_ug_g"
)

# SEC-ICP-MS signal columns used only for ROI / overlay.
metal_signal_map <- tribble(
  ~metal, ~signal_col,
  "Cd",   "Cd113",
  "Cu",   "Cu65",
  "Zn",   "Zn68",
  "Ni",   "Ni61",
  "Pb",   "Pb206"
)

target_tissues <- c("gill", "viscera", "muscle")
target_tissue_aliases <- c("gill", "gills", "vis", "viscera", "muscle")

# ROI-weighted Elastic Net settings.
# Soft guidance is recommended. Avoid penalty.factor = 0 for n is small.
use_cd_roi_weighted_models <- TRUE
roi_penalty_core       <- 0.25
roi_penalty_flank      <- 0.50
roi_penalty_background <- 1.00
roi_flank_width_min    <- 1.00

# Automatic ROI detection from training-set Cd-specific SEC-ICP-MS profiles.
roi_threshold_frac <- 0.20
roi_min_width_min  <- 0.40
roi_merge_gap_min  <- 0.40

# ============================================================
# 1. Helper functions
# ============================================================

clean_sample_id <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_squish() %>%
    stringr::str_replace_all("_", "-")
}

clean_tissue <- function(x) {
  x <- as.character(x) %>%
    stringr::str_squish() %>%
    stringr::str_to_lower()
  
  dplyr::case_when(
    x %in% c("gill", "gills") ~ "gill",
    x %in% c("vis", "viscera") ~ "viscera",
    x == "muscle" ~ "muscle",
    TRUE ~ x
  )
}

make_positive_pseudo <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) stop("No positive values available for pseudo-count.")
  min(x, na.rm = TRUE) / 2
}

get_dt_min <- function(time_min) {
  time_min <- sort(unique(time_min[is.finite(time_min)]))
  dt <- median(diff(time_min), na.rm = TRUE)
  if (!is.finite(dt) || dt <= 0) {
    stop("Cannot estimate a valid retention-time interval.")
  }
  dt
}

safe_r2 <- function(obs, pred) {
  cc <- complete.cases(obs, pred)
  obs <- obs[cc]
  pred <- pred[cc]
  
  if (length(obs) < 3) return(NA_real_)
  if (sd(obs) == 0 || sd(pred) == 0) return(NA_real_)
  
  cor(obs, pred)^2
}

metric_tbl <- function(obs, pred, model_name, feature_set, split_name) {
  tibble(
    model = model_name,
    feature_set = feature_set,
    split = split_name,
    RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
    MAE = mean(abs(obs - pred), na.rm = TRUE),
    R2 = safe_r2(obs, pred),
    n_test = sum(complete.cases(obs, pred))
  )
}

calc_metrics_simple <- function(obs, pred) {
  tibble(
    RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
    MAE = mean(abs(obs - pred), na.rm = TRUE),
    R2 = safe_r2(obs, pred),
    n = sum(complete.cases(obs, pred))
  )
}

drop_zero_variance_terms <- function(dat, x_terms) {
  x_terms <- intersect(x_terms, names(dat))
  x_terms <- x_terms[sapply(dat[x_terms], is.numeric)]
  
  if (length(x_terms) == 0) return(character(0))
  
  sds <- sapply(dat[x_terms], sd, na.rm = TRUE)
  x_terms[is.finite(sds) & sds > 0]
}

theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold")
    )
}

label_no_trailing0 <- function(x) {
  sapply(
    x,
    function(z) {
      out <- format(z, scientific = FALSE, trim = TRUE, digits = 10)
      if (grepl("\\.", out)) {
        out <- sub("0+$", "", out)
        out <- sub("\\.$", "", out)
      }
      out
    }
  )
}

make_axis_limits <- function(x, y) {
  z <- c(x, y)
  z <- z[is.finite(z) & z > 0]
  if (length(z) == 0) stop("No positive values for axis limits.")
  
  low <- 10^(floor(log10(min(z))) - 0.05)
  high <- 10^(ceiling(log10(max(z))) + 0.05)
  
  c(low, high)
}

smooth_vec <- function(x, k = 5) {
  x <- as.numeric(x)
  if (length(x) < k) return(x)
  y <- as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
  y[!is.finite(y)] <- x[!is.finite(y)]
  y
}

extract_rt_interval <- function(x_terms) {
  tibble(
    predictor = x_terms,
    rt_start = as.numeric(stringr::str_match(x_terms, "rt_(\\d+(?:\\.\\d+)?)_\\d+(?:\\.\\d+)?")[, 2]),
    rt_end   = as.numeric(stringr::str_match(x_terms, "rt_\\d+(?:\\.\\d+)?_(\\d+(?:\\.\\d+)?)")[, 2])
  ) %>%
    mutate(
      rt_mid = (rt_start + rt_end) / 2
    )
}


# ============================================================
# 2. Read metal concentration table
# ============================================================
metal_raw <- metal_raw %>%
  mutate(
    sample_id = clean_sample_id(sample),
    Species = as.character(Species),
    Tissue = clean_tissue(Tissue),
    Fraction = as.character(Fraction),
    rep = as.character(rep),
    ww_g = as.numeric(ww_g)
  )

missing_total_input_cols <- setdiff(
  c(metal_total_map$back_col, metal_total_map$accumulated_col),
  names(metal_raw)
)

if (length(missing_total_input_cols) > 0) {
  stop(
    "Missing metal input column(s) in metal file: ",
    paste(missing_total_input_cols, collapse = ", ")
  )
}

for (i in seq_len(nrow(metal_total_map))) {
  total_col <- metal_total_map$total_col[i]
  back_col <- metal_total_map$back_col[i]
  accumulated_col <- metal_total_map$accumulated_col[i]
  
  metal_raw[[total_col]] <-
    as.numeric(metal_raw[[back_col]]) +
    as.numeric(metal_raw[[accumulated_col]])
}

build_metal_target_data <- function(met, target_col) {
  if (!target_col %in% names(metal_raw)) {
    stop("Target concentration column not found: ", target_col)
  }
  
  out <- metal_raw %>%
    filter(stringr::str_to_lower(Fraction) == "total") %>%
    transmute(
      sample_id,
      Species,
      Tissue,
      rep,
      ww_g,
      metal = met,
      metal_conc = as.numeric(.data[[target_col]])
    ) %>%
    filter(
      Tissue %in% target_tissues,
      is.finite(metal_conc)
    ) %>%
    mutate(
      Tissue = factor(Tissue, levels = target_tissues),
      Species = factor(Species)
    ) %>%
    distinct()
  
  if (!any(is.finite(out$metal_conc) & out$metal_conc > 0)) {
    stop("No positive concentration values found for ", met)
  }
  
  pseudo <- make_positive_pseudo(out$metal_conc)
  
  out %>%
    mutate(
      metal_pseudo = pseudo,
      metal_log10 = log10(metal_conc + metal_pseudo)
    )
}

cd_target_data <- build_metal_target_data(
  met = target_metal,
  target_col = target_col
) %>%
  rename(
    Cd = metal_conc,
    Cd_pseudo = metal_pseudo,
    Cd_log10 = metal_log10
  )

cat("\nCd target data from metal concentration table:\n")
print(cd_target_data %>% count(Species, Tissue))

cat("\nCd concentration summary, unit = ug g^-1 wet weight:\n")
print(summary(cd_target_data$Cd))


# ============================================================
# 3. Read HPLC-UV data
# ============================================================

names(protein_raw)[1] <- "time_min"

protein_raw <- protein_raw %>%
  mutate(time_min = as.numeric(time_min)) %>%
  filter(is.finite(time_min))

uv_dt <- get_dt_min(protein_raw$time_min)

protein_long <- protein_raw %>%
  pivot_longer(
    cols = -time_min,
    names_to = "sample_id_raw",
    values_to = "uv"
  ) %>%
  mutate(
    sample_id = clean_sample_id(sample_id_raw),
    uv = as.numeric(uv)
  ) %>%
  filter(is.finite(time_min))

matched_ids <- intersect(unique(protein_long$sample_id), unique(cd_target_data$sample_id))

if (length(matched_ids) == 0) {
  stop("No matched sample IDs between HPLC-UV and metal concentration table.")
}

protein_long <- protein_long %>%
  filter(sample_id %in% matched_ids)

model_meta <- cd_target_data %>%
  filter(sample_id %in% matched_ids) %>%
  mutate(
    Tissue = factor(as.character(Tissue), levels = target_tissues),
    Species = droplevels(Species)
  )

cat("\nMatched UV-metal samples used for Cd modelling:\n")
print(model_meta %>% count(Species, Tissue))


# ============================================================
# 4. Read SEC-ICP-MS data for ROI prior only
# ============================================================

required_icpms_cols <- c("sample", "Species", "Tissue", "rep", "Time_sec",
                         target_icpms_signal)

missing_icpms_cols <- setdiff(required_icpms_cols, names(icpms_raw))

if (length(missing_icpms_cols) > 0) {
  stop(
    "Missing required column(s) in ICP-MS file: ",
    paste(missing_icpms_cols, collapse = ", ")
  )
}

icpms_long <- icpms_raw %>%
  mutate(
    sample_id = clean_sample_id(sample),
    Species = as.character(Species),
    Tissue = clean_tissue(Tissue),
    rep = as.character(rep),
    time_min = as.numeric(Time_sec) / 60
  ) %>%
  filter(
    Tissue %in% target_tissues,
    is.finite(time_min),
    time_min >= 0,
    time_min <= max_time
  )

for (cc in intersect(metal_signal_map$signal_col, names(icpms_long))) {
  icpms_long[[cc]] <- as.numeric(icpms_long[[cc]])
}

cat("\nSEC-ICP-MS data loaded for ROI prior only. These signals are not predictors.\n")
print(icpms_long %>% distinct(sample_id, Species, Tissue) %>% count(Species, Tissue))


# ============================================================
# 5. UV feature engineering
# ============================================================

make_clr_from_wide <- function(wide_raw) {
  feature_cols <- setdiff(names(wide_raw), "sample_id")
  pseudo <- make_positive_pseudo(as.matrix(wide_raw[, feature_cols]))
  
  x <- as.matrix(wide_raw[, feature_cols]) + pseudo
  x_prop <- x / rowSums(x)
  x_clr <- log(x_prop) - rowMeans(log(x_prop))
  
  out <- bind_cols(
    wide_raw %>% select(sample_id),
    as_tibble(x_clr, .name_repair = "minimal")
  )
  
  names(out)[-1] <- paste0("clr_", feature_cols)
  out
}

make_fixed_auc_features <- function(
    protein_long,
    bin_width = 1,
    max_time = 60
) {
  
  n_bins <- max_time / bin_width
  
  bin_long <- protein_long %>%
    filter(
      time_min >= 0,
      time_min <= max_time
    ) %>%
    mutate(
      uv_pos = pmax(uv, 0),
      rt_bin = pmin(floor(time_min / bin_width), n_bins - 1),
      rt_start = rt_bin * bin_width,
      rt_end = rt_start + bin_width,
      rt_bin_name = paste0(
        "rt_",
        sprintf("%02d", rt_start),
        "_",
        sprintf("%02d", rt_end)
      )
    ) %>%
    group_by(sample_id, rt_bin_name) %>%
    summarise(
      bin_auc = sum(uv_pos, na.rm = TRUE) * uv_dt,
      .groups = "drop"
    )
  
  wide_raw <- bin_long %>%
    select(sample_id, rt_bin_name, bin_auc) %>%
    pivot_wider(
      names_from = rt_bin_name,
      values_from = bin_auc,
      values_fill = 0
    )
  
  make_clr_from_wide(wide_raw)
}

make_quantile_bin_features <- function(
    protein_long,
    n_bins = 30
) {
  
  time_vec <- protein_long$time_min
  time_vec <- time_vec[is.finite(time_vec)]
  
  breaks <- quantile(
    time_vec,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE,
    names = FALSE
  )
  
  breaks <- unique(breaks)
  
  if (length(breaks) < 4) {
    stop("Too few unique quantile breaks.")
  }
  
  bin_long <- protein_long %>%
    filter(
      time_min >= min(breaks),
      time_min <= max(breaks)
    ) %>%
    mutate(
      uv_pos = pmax(uv, 0),
      rt_bin = cut(
        time_min,
        breaks = breaks,
        include.lowest = TRUE,
        labels = FALSE
      ),
      rt_bin_name = paste0("qbin_", sprintf("%02d", rt_bin))
    ) %>%
    group_by(sample_id, rt_bin_name) %>%
    summarise(
      bin_auc = sum(uv_pos, na.rm = TRUE) * uv_dt,
      .groups = "drop"
    )
  
  wide_raw <- bin_long %>%
    select(sample_id, rt_bin_name, bin_auc) %>%
    pivot_wider(
      names_from = rt_bin_name,
      values_from = bin_auc,
      values_fill = 0
    )
  
  make_clr_from_wide(wide_raw)
}

denoise_wavelet_one_curve <- function(y) {
  y <- as.numeric(y)
  y[!is.finite(y)] <- 0
  
  if (!has_waveslim) {
    return(pmax(y, 0))
  }
  
  n <- length(y)
  n_pad <- 2^ceiling(log2(n))
  
  y_reflect <- c(y, rev(y))
  y_pad <- rep(y_reflect, length.out = n_pad)
  
  n_level <- min(4, floor(log2(n_pad)) - 1)
  
  if (n_level < 1) {
    return(pmax(y, 0))
  }
  
  wt <- waveslim::dwt(
    y_pad,
    wf = "la8",
    n.levels = n_level
  )
  
  detail_names <- grep("^d", names(wt), value = TRUE)
  
  for (nm in detail_names) {
    sigma <- median(abs(wt[[nm]]), na.rm = TRUE) / 0.6745
    
    if (!is.finite(sigma) || sigma <= 0) {
      next
    }
    
    threshold <- sigma * sqrt(2 * log(length(wt[[nm]])))
    wt[[nm]] <- sign(wt[[nm]]) * pmax(abs(wt[[nm]]) - threshold, 0)
  }
  
  y_hat <- waveslim::idwt(wt)[seq_len(n)]
  
  pmax(y_hat, 0)
}

make_wavelet_auc_features <- function(
    protein_long,
    bin_width = 1,
    max_time = 60
) {
  
  if (!has_waveslim) {
    warning("Package waveslim is not installed. wavelet_clr will equal non-denoised fixed_auc_clr.")
  }
  
  denoised_long <- protein_long %>%
    arrange(sample_id, time_min) %>%
    group_by(sample_id) %>%
    mutate(
      uv = denoise_wavelet_one_curve(uv)
    ) %>%
    ungroup()
  
  make_fixed_auc_features(
    protein_long = denoised_long,
    bin_width = bin_width,
    max_time = max_time
  )
}

make_peak_features <- function(
    protein_long,
    minpeakdistance = 8,
    threshold_frac = 0.03,
    window_min = 0.35
) {
  
  if (!has_pracma) {
    warning("Package pracma is not installed. Skipping peak_clr features.")
    return(NULL)
  }
  
  mean_curve <- protein_long %>%
    group_by(time_min) %>%
    summarise(
      mean_uv = mean(pmax(uv, 0), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(time_min)
  
  peak_mat <- pracma::findpeaks(
    mean_curve$mean_uv,
    minpeakdistance = minpeakdistance,
    threshold = max(mean_curve$mean_uv, na.rm = TRUE) * threshold_frac
  )
  
  if (is.null(peak_mat) || nrow(peak_mat) < 2) {
    warning("Too few consensus peaks detected. Skipping peak_clr features.")
    return(NULL)
  }
  
  peak_index <- peak_mat[, 2]
  peak_time <- mean_curve$time_min[peak_index]
  
  peak_long <- map_dfr(seq_along(peak_time), function(i) {
    lo <- peak_time[i] - window_min
    hi <- peak_time[i] + window_min
    
    protein_long %>%
      filter(
        time_min >= lo,
        time_min <= hi
      ) %>%
      mutate(
        uv_pos = pmax(uv, 0)
      ) %>%
      group_by(sample_id) %>%
      summarise(
        peak_auc = sum(uv_pos, na.rm = TRUE) * uv_dt,
        .groups = "drop"
      ) %>%
      mutate(
        peak_name = paste0(
          "peak_",
          sprintf("%02d", i),
          "_",
          round(peak_time[i], 2),
          "min"
        )
      )
  })
  
  wide_raw <- peak_long %>%
    select(sample_id, peak_name, peak_auc) %>%
    pivot_wider(
      names_from = peak_name,
      values_from = peak_auc,
      values_fill = 0
    )
  
  make_clr_from_wide(wide_raw)
}


# ============================================================
# 6. SEC-ICP-MS-informed ROI utilities
# ============================================================

derive_roi_from_icpms <- function(
    icpms_long,
    training_sample_ids,
    signal_col = "Cd113",
    threshold_frac = 0.20,
    min_width_min = 0.40,
    merge_gap_min = 0.40,
    max_time = 60
) {
  
  if (!signal_col %in% names(icpms_long)) {
    stop("Signal column not found for ROI derivation: ", signal_col)
  }
  
  mean_curve <- icpms_long %>%
    filter(
      sample_id %in% training_sample_ids,
      Tissue %in% target_tissues,
      time_min >= 0,
      time_min <= max_time
    ) %>%
    mutate(
      signal = pmax(as.numeric(.data[[signal_col]]), 0)
    ) %>%
    group_by(time_min) %>%
    summarise(
      mean_signal = mean(signal, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(time_min) %>%
    mutate(
      smooth_signal = smooth_vec(mean_signal, k = 5)
    )
  
  if (nrow(mean_curve) == 0 || max(mean_curve$smooth_signal, na.rm = TRUE) <= 0) {
    warning("No valid ICP-MS signal available for ROI derivation. Using no ROI prior.")
    return(tibble(roi_name = character(), start_min = numeric(), end_min = numeric()))
  }
  
  threshold <- max(mean_curve$smooth_signal, na.rm = TRUE) * threshold_frac
  
  flag <- mean_curve$smooth_signal >= threshold
  
  if (!any(flag, na.rm = TRUE)) {
    warning("No ROI passed the threshold. Using no ROI prior.")
    return(tibble(roi_name = character(), start_min = numeric(), end_min = numeric()))
  }
  
  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  
  roi_tbl <- tibble(
    flag = r$values,
    idx_start = starts,
    idx_end = ends
  ) %>%
    filter(flag) %>%
    transmute(
      start_min = mean_curve$time_min[idx_start],
      end_min = mean_curve$time_min[idx_end]
    ) %>%
    mutate(
      width = end_min - start_min
    ) %>%
    filter(width >= min_width_min)
  
  if (nrow(roi_tbl) == 0) {
    warning("All detected ROIs were narrower than min_width_min. Using no ROI prior.")
    return(tibble(roi_name = character(), start_min = numeric(), end_min = numeric()))
  }
  
  roi_tbl <- roi_tbl %>%
    arrange(start_min)
  
  # Merge ROIs separated by small gaps.
  merged <- list()
  current <- roi_tbl[1, ]
  
  if (nrow(roi_tbl) > 1) {
    for (i in 2:nrow(roi_tbl)) {
      gap <- roi_tbl$start_min[i] - current$end_min
      if (is.finite(gap) && gap <= merge_gap_min) {
        current$end_min <- roi_tbl$end_min[i]
        current$width <- current$end_min - current$start_min
      } else {
        merged[[length(merged) + 1]] <- current
        current <- roi_tbl[i, ]
      }
    }
  }
  
  merged[[length(merged) + 1]] <- current
  
  bind_rows(merged) %>%
    mutate(
      roi_name = paste0(
        "Cd_ROI_",
        sprintf("%02d", row_number()),
        "_",
        round(start_min, 2),
        "_",
        round(end_min, 2),
        "min"
      )
    ) %>%
    select(roi_name, start_min, end_min)
}

make_cd_roi_penalty_factor <- function(
    x_terms,
    roi_tbl,
    core_penalty = roi_penalty_core,
    flank_penalty = roi_penalty_flank,
    background_penalty = roi_penalty_background,
    flank_width_min = roi_flank_width_min
) {
  
  pf <- rep(background_penalty, length(x_terms))
  names(pf) <- x_terms
  
  if (is.null(roi_tbl) || nrow(roi_tbl) == 0) {
    return(pf)
  }
  
  rt_tbl <- extract_rt_interval(x_terms)
  
  for (i in seq_len(nrow(roi_tbl))) {
    roi_start <- as.numeric(roi_tbl$start_min[i])
    roi_end <- as.numeric(roi_tbl$end_min[i])
    
    core_idx <- which(
      is.finite(rt_tbl$rt_start) &
        is.finite(rt_tbl$rt_end) &
        rt_tbl$rt_start < roi_end &
        rt_tbl$rt_end > roi_start
    )
    
    flank_start <- roi_start - flank_width_min
    flank_end <- roi_end + flank_width_min
    
    flank_idx <- which(
      is.finite(rt_tbl$rt_start) &
        is.finite(rt_tbl$rt_end) &
        rt_tbl$rt_start < flank_end &
        rt_tbl$rt_end > flank_start
    )
    
    pf[flank_idx] <- pmin(pf[flank_idx], flank_penalty)
    pf[core_idx] <- pmin(pf[core_idx], core_penalty)
  }
  
  pf
}


# ============================================================
# 7. Build feature sets
# ============================================================

feature_sets <- list(
  fixed_auc_clr = make_fixed_auc_features(
    protein_long = protein_long,
    bin_width = bin_width,
    max_time = max_time
  ),
  
  quantile_clr = make_quantile_bin_features(
    protein_long = protein_long,
    n_bins = n_quantile_bins
  ),
  
  wavelet_clr = make_wavelet_auc_features(
    protein_long = protein_long,
    bin_width = bin_width,
    max_time = max_time
  )
)

peak_features <- make_peak_features(protein_long)

if (!is.null(peak_features)) {
  feature_sets$peak_clr <- peak_features
}

# ROI-weighted feature sets use the same UV predictors but different glmnet penalties.
if (use_cd_roi_weighted_models) {
  feature_sets$fixed_auc_clr_roi_weighted <- feature_sets$fixed_auc_clr
  feature_sets$wavelet_clr_roi_weighted <- feature_sets$wavelet_clr
}

cat("\nUV feature sets created:\n")
print(names(feature_sets))


# ============================================================
# 8. Modelling functions
# ============================================================

prepare_xy <- function(feature_df, target_meta) {
  if (!is.data.frame(feature_df)) {
    stop("feature_df must be a data.frame/tibble, but got: ", class(feature_df)[1])
  }
  
  dat <- target_meta %>%
    inner_join(feature_df, by = "sample_id") %>%
    filter(is.finite(Cd_log10))
  
  meta_cols <- c(
    "sample_id",
    "Species",
    "Tissue",
    "rep",
    "ww_g",
    "metal",
    "Cd",
    "Cd_pseudo",
    "Cd_log10"
  )
  
  x_terms <- setdiff(names(dat), meta_cols)
  x_terms <- drop_zero_variance_terms(dat, x_terms)
  
  if (length(x_terms) < 1) {
    stop("No usable predictors remain after variance filtering.")
  }
  
  dat <- dat %>%
    filter(
      if_all(all_of(x_terms), ~ is.finite(.x))
    )
  
  list(
    dat = dat,
    x_terms = x_terms
  )
}

fit_predict_mean <- function(train_dat, test_dat) {
  pred <- rep(mean(train_dat$Cd_log10, na.rm = TRUE), nrow(test_dat))
  
  list(
    pred = pred,
    coef_tbl = tibble(
      predictor = character(0),
      beta = numeric(0),
      alpha = numeric(0),
      lambda = numeric(0),
      penalty_factor = numeric(0)
    ),
    alpha = NA_real_,
    lambda = NA_real_
  )
}

fit_predict_lm_single <- function(train_dat, test_dat, x_term) {
  train_lm <- train_dat %>%
    transmute(
      Cd_log10 = Cd_log10,
      .x = .data[[x_term]]
    )
  
  test_lm <- test_dat %>%
    transmute(
      .x = .data[[x_term]]
    )
  
  fit <- lm(Cd_log10 ~ .x, data = train_lm)
  
  pred <- as.numeric(
    predict(
      fit,
      newdata = test_lm
    )
  )
  
  beta <- coef(fit)
  
  coef_tbl <- tibble(
    predictor = x_term,
    beta = as.numeric(beta[".x"]),
    alpha = NA_real_,
    lambda = NA_real_,
    penalty_factor = NA_real_
  )
  
  list(
    pred = pred,
    coef_tbl = coef_tbl,
    alpha = NA_real_,
    lambda = NA_real_
  )
}

fit_predict_glmnet <- function(
    train_dat,
    test_dat,
    x_terms,
    penalty_factor = NULL,
    alpha_grid = c(0.25, 0.5, 0.75, 1),
    nfolds_inner = 5
) {
  
  if (nrow(train_dat) < 4 || length(unique(train_dat$Cd_log10)) < 3) {
    warning("Training set too small or y has too little variation. Using mean predictor.")
    return(fit_predict_mean(train_dat, test_dat))
  }
  
  x_train <- as.matrix(train_dat[, x_terms, drop = FALSE])
  x_test  <- as.matrix(test_dat[, x_terms, drop = FALSE])
  y_train <- train_dat$Cd_log10
  
  sds <- apply(x_train, 2, sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  
  x_train <- x_train[, keep, drop = FALSE]
  x_test  <- x_test[, keep, drop = FALSE]
  kept_terms <- colnames(x_train)
  
  if (ncol(x_train) < 1) {
    warning("No predictors left after train-only variance filtering. Using mean predictor.")
    return(fit_predict_mean(train_dat, test_dat))
  }
  
  if (ncol(x_train) == 1) {
    message("Only one predictor after filtering; using lm() instead of glmnet().")
    return(
      fit_predict_lm_single(
        train_dat = train_dat,
        test_dat = test_dat,
        x_term = kept_terms[1]
      )
    )
  }
  
  if (is.null(penalty_factor)) {
    pf_kept <- rep(1, ncol(x_train))
    names(pf_kept) <- kept_terms
  } else {
    pf_kept <- penalty_factor[kept_terms]
    pf_kept[!is.finite(pf_kept) | is.na(pf_kept)] <- 1
  }
  
  inner_nfolds <- min(nfolds_inner, nrow(train_dat))
  
  if (inner_nfolds < 3) {
    warning("Too few training samples for inner CV. Using mean predictor.")
    return(fit_predict_mean(train_dat, test_dat))
  }
  
  cvfits <- map(alpha_grid, function(a) {
    cv.glmnet(
      x = x_train,
      y = y_train,
      alpha = a,
      family = "gaussian",
      standardize = TRUE,
      nfolds = inner_nfolds,
      penalty.factor = pf_kept
    )
  })
  
  cv_min <- map_dbl(cvfits, ~ min(.x$cvm, na.rm = TRUE))
  best_i <- which.min(cv_min)
  
  pred <- as.numeric(
    predict(
      cvfits[[best_i]],
      newx = x_test,
      s = "lambda.1se"
    )
  )
  
  coef_mat <- coef(cvfits[[best_i]], s = "lambda.1se")
  
  coef_tbl <- tibble(
    predictor = rownames(coef_mat),
    beta = as.numeric(coef_mat[, 1]),
    alpha = alpha_grid[best_i],
    lambda = cvfits[[best_i]]$lambda.1se
  ) %>%
    filter(
      predictor != "(Intercept)",
      beta != 0
    ) %>%
    mutate(
      penalty_factor = pf_kept[predictor]
    )
  
  list(
    pred = pred,
    coef_tbl = coef_tbl,
    alpha = alpha_grid[best_i],
    lambda = cvfits[[best_i]]$lambda.1se
  )
}


# ============================================================
# 9. Leave-One-Species-Out CV for Cd
# ============================================================

run_loso_one_feature_set <- function(feature_df, feature_name) {
  px <- prepare_xy(
    feature_df = feature_df,
    target_meta = model_meta
  )
  
  dat <- px$dat
  x_terms <- px$x_terms
  
  species_levels <- levels(droplevels(dat$Species))
  
  if (length(species_levels) < 2) {
    stop("LOSO requires at least two species.")
  }
  
  all_pred <- list()
  all_metric <- list()
  all_coef <- list()
  all_roi <- list()
  
  for (sp in species_levels) {
    train_dat <- dat %>%
      filter(Species != sp)
    
    test_dat <- dat %>%
      filter(Species == sp)
    
    split_name <- paste0(
      "Train other species -> test ",
      as.character(sp)
    )
    
    cat("\nRunning ", feature_name, " | ", split_name, "\n", sep = "")
    
    penalty_factor <- NULL
    roi_tbl_split <- tibble(
      roi_name = character(),
      start_min = numeric(),
      end_min = numeric()
    )
    
    if (stringr::str_detect(feature_name, "_roi_weighted$")) {
      roi_tbl_split <- derive_roi_from_icpms(
        icpms_long = icpms_long,
        training_sample_ids = train_dat$sample_id,
        signal_col = target_icpms_signal,
        threshold_frac = roi_threshold_frac,
        min_width_min = roi_min_width_min,
        merge_gap_min = roi_merge_gap_min,
        max_time = max_time
      )
      
      penalty_factor <- make_cd_roi_penalty_factor(
        x_terms = x_terms,
        roi_tbl = roi_tbl_split
      )
      
      all_roi[[length(all_roi) + 1]] <- roi_tbl_split %>%
        mutate(
          feature_set = feature_name,
          split = split_name,
          left_out_species = as.character(sp)
        )
      
      cat("Training-set Cd-binding ROI from SEC-ICP-MS:\n")
      print(roi_tbl_split)
    }
    
    en <- fit_predict_glmnet(
      train_dat = train_dat,
      test_dat = test_dat,
      x_terms = x_terms,
      penalty_factor = penalty_factor
    )
    
    pred_en <- test_dat %>%
      transmute(
        sample_id,
        Species,
        Tissue,
        rep,
        Cd,
        Cd_log10,
        pred_log10 = en$pred,
        model = "Elastic Net",
        feature_set = feature_name,
        split = split_name,
        alpha = en$alpha,
        lambda = en$lambda
      )
    
    all_pred[[length(all_pred) + 1]] <- pred_en
    
    all_metric[[length(all_metric) + 1]] <- metric_tbl(
      obs = test_dat$Cd_log10,
      pred = en$pred,
      model_name = "Elastic Net",
      feature_set = feature_name,
      split_name = split_name
    ) %>%
      mutate(
        alpha = en$alpha,
        lambda = en$lambda
      )
    
    all_coef[[length(all_coef) + 1]] <- en$coef_tbl %>%
      mutate(
        feature_set = feature_name,
        split = split_name
      )
  }
  
  list(
    pred = bind_rows(all_pred),
    metrics_by_split = bind_rows(all_metric),
    coef = bind_rows(all_coef),
    roi = bind_rows(all_roi)
  )
}

loso_results <- imap(feature_sets, run_loso_one_feature_set)

loso_pred <- loso_results %>%
  map("pred") %>%
  bind_rows()

loso_metrics_by_split <- loso_results %>%
  map("metrics_by_split") %>%
  bind_rows() %>%
  arrange(feature_set, model, split)

elastic_net_coef <- loso_results %>%
  map("coef") %>%
  bind_rows()

roi_by_split <- loso_results %>%
  map("roi") %>%
  bind_rows()

loso_metrics_overall <- loso_pred %>%
  group_by(feature_set, model) %>%
  summarise(
    RMSE = sqrt(mean((Cd_log10 - pred_log10)^2, na.rm = TRUE)),
    MAE = mean(abs(Cd_log10 - pred_log10), na.rm = TRUE),
    R2 = safe_r2(Cd_log10, pred_log10),
    n = sum(complete.cases(Cd_log10, pred_log10)),
    .groups = "drop"
  ) %>%
  arrange(RMSE)

cat("\nLOSO metrics by split:\n")
print(loso_metrics_by_split)

cat("\nOverall LOSO metrics:\n")
print(loso_metrics_overall)

best_pair <- loso_metrics_overall %>%
  arrange(RMSE) %>%
  slice(1) %>%
  select(feature_set, model)

loso_pred_plot <- loso_pred %>%
  mutate(
    pred_Cd = 10^pred_log10 - unique(model_meta$Cd_pseudo),
    pred_Cd = pmax(pred_Cd, 0)
  )

# ============================================================
# ============================================================
# 10. Cd prediction plots and selected features
# ============================================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ============================================================
# 10.3 Prepare plotting data for best Cd model
# ============================================================

plot_df_cd_best <- loso_pred_plot %>%
  semi_join(
    best_pair,
    by = c("feature_set", "model")
  ) %>%
  mutate(
    Tissue = factor(
      as.character(Tissue),
      levels = c("gill", "viscera", "muscle")
    ),
    Species = factor(
      as.character(Species),
      levels = c("A. kagoshimensis", "T. granosa")
    )
  )

if (nrow(plot_df_cd_best) == 0) {
  stop("No Cd prediction data available for the best model.")
}

axis_lims_cd <- make_axis_limits(
  plot_df_cd_best$Cd,
  plot_df_cd_best$pred_Cd
)

ref_line_cd <- tibble(
  x = 10^seq(
    log10(axis_lims_cd[1]),
    log10(axis_lims_cd[2]),
    length.out = 200
  ),
  y_1to1 = x,
  y_2x = 2 * x,
  y_halfx = 0.5 * x
)

log_breaks_cd <- c(
  0.01, 
  0.1,
  1,
  10
)

log_breaks_cd_use <- log_breaks_cd[
  log_breaks_cd >= axis_lims_cd[1] &
    log_breaks_cd <= axis_lims_cd[2]
]

# ------------------------------------------------------------
# Aesthetic settings
# Species = color
# Tissue  = shape
# ------------------------------------------------------------


shape_vals_tissue <- c(
  "gill"    = 21,  # circle, fillable
  "viscera" = 24,  # triangle, fillable
  "muscle"  = 22   # square, fillable
)

species_cols <- c(
  "A. kagoshimensis" = "#F2A3A8",
  "T. granosa"       = "#4F7FD9"
)

# ------------------------------------------------------------
# Statistics label for the best Cd model
# ------------------------------------------------------------

best_metric_cd <- loso_metrics_overall %>%
  semi_join(
    best_pair,
    by = c("feature_set", "model")
  ) %>%
  slice(1)

stat_label_cd <- tibble(
  label_x = 10^(
    log10(axis_lims_cd[1]) +
      0.06 * (log10(axis_lims_cd[2]) - log10(axis_lims_cd[1]))
  ),
  label_y = 10^(
    log10(axis_lims_cd[2]) -
      0.08 * (log10(axis_lims_cd[2]) - log10(axis_lims_cd[1]))
  ),
  label = paste0(
    "RMSE = ", sprintf("%.2f", best_metric_cd$RMSE),
    "\nMAE = ", sprintf("%.2f", best_metric_cd$MAE),
    "\nR² = ", ifelse(is.na(best_metric_cd$R2), "NA", sprintf("%.2f", best_metric_cd$R2))
  )
)

# ------------------------------------------------------------
# Species text labels inside plot
# ------------------------------------------------------------
species_label_pos <- tribble(
  ~Species,             ~x_prop, ~y_prop,
  "A. kagoshimensis",   0.5,    0.82,
  "T. granosa",         0.67,    0.66
)

species_label_cd <- species_label_pos %>%
  mutate(
    Species = factor(
      Species,
      levels = c("A. kagoshimensis", "T. granosa")
    ),
     label_x = 10^(
      log10(axis_lims_cd[1]) +
        x_prop * (log10(axis_lims_cd[2]) - log10(axis_lims_cd[1]))
    ),
    label_y = 10^(
      log10(axis_lims_cd[1]) +
        y_prop * (log10(axis_lims_cd[2]) - log10(axis_lims_cd[1]))
    ),
    label = case_when(
      Species == "A. kagoshimensis" ~ "italic('A. kagoshimensis')",
      Species == "T. granosa"       ~ "italic('T. granosa')"
    )
  )

# ============================================================
# 10.4 Plot Cd observed vs predicted
# ============================================================
p_cd_loso <- ggplot(
  plot_df_cd_best,
  aes(
    x = Cd,
    y = pred_Cd
  )
) +
  geom_line(
    data = ref_line_cd,
    aes(x = x, y = y_1to1),
    inherit.aes = FALSE,
    linewidth = 0.7,
    color = "black"
  ) +
  geom_line(
    data = ref_line_cd,
    aes(x = x, y = y_2x),
    inherit.aes = FALSE,
    linetype = 2,
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_line(
    data = ref_line_cd,
    aes(x = x, y = y_halfx),
    inherit.aes = FALSE,
    linetype = 2,
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_point(
    aes(
      fill = Species,
      shape = Tissue
    ),
    color = "black",
    size = 3.2,
    stroke = 0.7,
    alpha = 0.95
  ) +
  geom_text(
    data = stat_label_cd,
    aes(
      x = label_x,
      y = label_y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.2
  ) +
  geom_text(
    data = species_label_cd,
    aes(
      x = label_x,
      y = label_y,
      label = label,
      color = Species
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.2,
    parse = TRUE,
    show.legend = FALSE
  )+
  scale_x_log10(
    limits = axis_lims_cd,
    breaks = log_breaks_cd_use,
    labels = label_no_trailing0
  ) +
  scale_y_log10(
    limits = axis_lims_cd,
    breaks = log_breaks_cd_use,
    labels = label_no_trailing0
  ) +
  scale_fill_manual(
    values = species_cols,
    guide = "none"
  ) +
  scale_color_manual(
    values = species_cols,
    guide = "none"
  ) +
  scale_shape_manual(
    values = shape_vals_tissue,
    name = "Tissue",
    breaks = c("gill", "viscera", "muscle")
  )+
  theme_pub(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.82, 0.23),
    legend.justification = c(0.5, 0.5),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )+
  labs(
    x = expression("Observed Cd conc. ("*mu*"g g"^-1*" ww)"),
    y = expression("Predicted Cd conc. ("*mu*"g g"^-1*" ww)")
  )

print(p_cd_loso)


# ============================================================
# 10.5 Prepare Elastic Net selected-feature table
# ============================================================

selected_en_features <- elastic_net_coef %>%
  filter(beta != 0) %>%
  group_by(
    feature_set,
    predictor
  ) %>%
  summarise(
    selected_n = n(),
    beta_mean = mean(beta, na.rm = TRUE),
    beta_abs_mean = mean(abs(beta), na.rm = TRUE),
    penalty_factor_mean = mean(penalty_factor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    desc(selected_n),
    desc(beta_abs_mean)
  )

readr::write_csv(
  selected_en_features,
  "Cd_ElasticNet_selected_features_summary.csv"
)

selected_en_rt <- selected_en_features %>%
  mutate(
    rt_start = as.numeric(
      stringr::str_match(
        predictor,
        "rt_(\\d+)_\\d+"
      )[, 2]
    ),
    rt_end = as.numeric(
      stringr::str_match(
        predictor,
        "rt_\\d+_(\\d+)"
      )[, 2]
    ),
    rt_mid = (rt_start + rt_end) / 2
  ) %>%
  filter(
    is.finite(rt_mid)
  )

best_feature <- best_pair$feature_set[1]

selected_en_rt_best <- selected_en_rt %>%
  filter(
    feature_set == best_feature
  )

if (nrow(selected_en_rt_best) == 0) {
  warning("No retention-time-resolved Elastic Net coefficients for the best feature set.")
}


# ============================================================
# 10.6 Overlay of selected Elastic Net coefficients and
#      tissue-specific SEC-ICP-MS Cd profiles
# ============================================================
sec_signal_col <- "Cd113"  # measured total Cd113 signal

coef_max <- max(selected_en_rt_best$beta_abs_mean, na.rm = TRUE)
sec_scale <- coef_max * 0.85

sec_overlay_df <- icpms_long %>%
  filter(
    Tissue %in% c("gill", "viscera", "muscle"),
    time_min >= 0,
    time_min <= max_time
  ) %>%
  mutate(
    sec_signal = pmax(as.numeric(.data[[sec_signal_col]]), 0),
    Tissue_plot = factor(
      Tissue,
      levels = c("gill", "viscera", "muscle"),
      labels = c("gill", "viscera", "muscle")
    )
  ) %>%
  group_by(Tissue_plot, time_min) %>%
  summarise(
    sec_mean = mean(sec_signal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Tissue_plot) %>%
  mutate(
    sec_norm = sec_mean / max(sec_mean, na.rm = TRUE),
    sec_plot = sec_norm * sec_scale
  ) %>%
  ungroup()

p_selected_rt_overlay <- ggplot() +
  geom_line(
    data = sec_overlay_df,
    aes(
      x = time_min,
      y = sec_plot,
      color = Tissue_plot
    ),
    linewidth = 0.4,
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_col(
    data = selected_en_rt_best,
    aes(
      x = rt_mid,
      y = beta_abs_mean
    ),
    width = 0.75,
    fill =  "#AB5B64",
    color = "white",
    linewidth = 0.15,
    alpha = 0.95
  ) +
  scale_color_manual(
    values = c(
      "gill" = "#F2A3A8",
      "viscera" = "#8FA8E8",
      "muscle" = "#3F5D68"
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, max_time, 10),
    limits = c(0, max_time),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, coef_max * 1.08),
    expand = expansion(mult = c(0, 0.04))
  ) +
  theme_pub(base_size = 12) +
  theme(
    legend.position = "none"
  ) +
  labs(
    x = "Retention time (min)",
    y = "Absolute Elastic Net coefficient"
  )

print(p_selected_rt_overlay)


# ============================================================
# 10.7 Optional: combined Cd figure
#     Keep this separate so individual panels remain easy to adjust.
# ============================================================
p1 <- p_cd_loso +
  theme(
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    axis.title.x = element_text(margin = margin(t = 8))
  )

p2 <- p_selected_rt_overlay +
  theme(
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    axis.title.x = element_text(margin = margin(t = 8))
  )

aligned_plots <- cowplot::align_plots(
  p1, p2,
  align = "h",
  axis = "b"
)

p_cd_combined <- cowplot::plot_grid(
  aligned_plots[[1]],
  aligned_plots[[2]],
  ncol = 2,
  labels = c("a", "b"),
  label_size = 16,
  label_fontface = "bold",
  rel_widths = c(1, 1)
)

print(p_cd_combined)
  
ggsave(
  filename = "Cd_mode2-1.png",
  width = 564/90,
  height =327/90,
  dpi = 900
)

ggsave(
  filename = "Cd_mode2-1.pdf",
  width = 564/90,
  height =327/90
)
