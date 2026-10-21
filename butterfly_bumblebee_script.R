suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readxl)
  library(stringr)
  library(lubridate)
  library(janitor)
  library(bombLepiSurv)
  library(DBI)
  library(glmmTMB)
  library(gllvm)
})
set.seed(123)

input_dir <- Sys.getenv("ANALYSIS_INPUT_DIR", unset = "data")
output_dir <- Sys.getenv("ANALYSIS_OUTPUT_DIR", unset = "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
traits_path <- file.path(input_dir, "Traits compilation 2026 bumblebees butterflies.xlsx")
landscape_path <- file.path(input_dir, "spatial_dat.rds")

abund_locality_col <- "Flate"
clim_locality_col <- "locality"
abund_meta_cols <- c("Flate", "Year", "Region", "Habitattype", "Blomsterdekke")
clim_date_col <- "date"
clim_temp_col <- "locality_tg"
clim_prec_col <- "locality_rr"
year_min <- 2011L
year_max <- 2025L
anomaly_reference_years <- NULL
min_days_per_window <- 20L
min_species_occurrences <- 10L

required_file <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  path
}
z_score <- function(x) {
  sx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / sx)
}


# ---- Shared functions ------------------------------------------------------
get_species_cols <- function(dat, meta_cols) {
  dat %>%
    select(where(is.numeric)) %>%
    select(-any_of(meta_cols)) %>%
    names()
}

prepare_abundance <- function(dat, taxon_group, meta_cols) {
  dat <- dat %>%
    as_tibble() %>%
    clean_names()

  meta_clean <- janitor::make_clean_names(meta_cols)
  locality_col_clean <- janitor::make_clean_names(abund_locality_col)
  species_cols <- get_species_cols(dat, meta_clean)

  if (!locality_col_clean %in% names(dat)) {
    stop("Locality column '", locality_col_clean, "' not found in ", taxon_group, " data.")
  }

  if (length(species_cols) == 0) {
    stop("No species columns found for ", taxon_group, ". Check abund_meta_cols.")
  }

  dat %>%
    mutate(
      taxon_group = taxon_group,
      total_abundance = rowSums(across(all_of(species_cols)), na.rm = TRUE)
    ) %>%
    rename(
      locality_id_abund = all_of(locality_col_clean),
      year = year,
      region = region,
      habitat_type = habitattype,
      flower_cover = blomsterdekke
    ) %>%
    mutate(
      taxon_group = factor(taxon_group, levels = c("Butterfly", "Bumblebee")),
      locality_id_abund = as.character(locality_id_abund),
      year = as.integer(year),
      region = factor(region),
      habitat_type = factor(habitat_type),
      flower_cover = as.numeric(flower_cover)
    )
}

make_climate_window <- function(
  clim,
  window_name,
  start_mmdd,
  end_mmdd,
  offset_year = 0,
  hot_threshold = 20
) {
  start_month <- as.integer(substr(start_mmdd, 1, 2))
  start_day   <- as.integer(substr(start_mmdd, 4, 5))
  end_month   <- as.integer(substr(end_mmdd, 1, 2))
  end_day     <- as.integer(substr(end_mmdd, 4, 5))

  clim_window <- clim %>%
    ungroup() %>%
    mutate(
      climate_year = lubridate::year(date),
      target_year = climate_year - offset_year,
      window_start = lubridate::make_date(
        year = climate_year,
        month = start_month,
        day = start_day
      ),
      window_end = lubridate::make_date(
        year = climate_year,
        month = end_month,
        day = end_day
      )
    ) %>%
    filter(date >= window_start, date <= window_end)

  clim_window %>%
    group_by(locality_id_clim, target_year) %>%
    summarise(
      first_date_used = min(date, na.rm = TRUE),
      last_date_used = max(date, na.rm = TRUE),
      n_rows = n(),
      n_dates = n_distinct(date),
      n_unique_temp = n_distinct(temp_mean),

      temp_mean_window = mean(temp_mean, na.rm = TRUE),
      temp_min_daily_mean = min(temp_mean, na.rm = TRUE),
      temp_max_daily_mean = max(temp_mean, na.rm = TRUE),
      temp_sd_daily_mean = sd(temp_mean, na.rm = TRUE),

      precip_sum = sum(precip, na.rm = TRUE),
      precip_mean_window = mean(precip, na.rm = TRUE),
      precip_max_daily = max(precip, na.rm = TRUE),

      gdd5 = sum(pmax(temp_mean - 5, 0), na.rm = TRUE),
      wet_days_1mm = sum(precip >= 1, na.rm = TRUE),
      dry_days_1mm = sum(precip < 1, na.rm = TRUE),

      hot_threshold_used = hot_threshold,
      hot_days = sum(temp_mean >= hot_threshold, na.rm = TRUE),

      frost_days_0c = sum(temp_mean < 0, na.rm = TRUE),

      .groups = "drop"
    ) %>%
    filter(n_dates >= min_days_per_window) %>%
    rename(year = target_year) %>%
    rename_with(
      .cols = -c(locality_id_clim, year),
      .fn = ~ paste(window_name, .x, sep = "__")
    )
}


add_climate_anomalies <- function(clim_summary, ref_years = NULL) { # grouping by loc, subtracting mean and z-scoring
  climate_vars <- clim_summary %>%
    select(where(is.numeric)) %>%
    select(-year) %>%
    names()

  if (!is.null(ref_years)) {
    ref_dat <- clim_summary %>%
      filter(year >= ref_years[1], year <= ref_years[2])
  } else {
    ref_dat <- clim_summary
  }

  ref_means <- ref_dat %>%
    group_by(locality_id_clim) %>%
    summarise(
      across(
        all_of(climate_vars),
        ~ mean(.x, na.rm = TRUE),
        .names = "{.col}__mean_ref"
      ),
      .groups = "drop"
    )

  out <- clim_summary %>%
    left_join(ref_means, by = "locality_id_clim")

  for (v in climate_vars) {
    ref_v <- paste0(v, "__mean_ref")
    anom_v <- paste0(v, "__anom")
    z_v <- paste0(v, "__anom__z")

    out[[anom_v]] <- out[[v]] - out[[ref_v]]
    out[[z_v]] <- z_score(out[[anom_v]]) # z-scoring climate vars
  }

  out
}



# ---- Load observations and climate -----------------------------------------

# A preconfigured, authenticated bombLepiSurv connection is required.
humlesommerfConnect()
fetch_taxon <- function(type) {
  grass <- getAllData(type = type, habitat = "Gressmark", year = year_min:year_max, language = "Latin")
  forest <- getAllData(type = type, habitat = "Skogsmark", year = year_min:year_max, language = "Latin")
  all <- bind_rows(grass, forest)
  metadata <- names(all)[seq_len(match("Periode", names(all)) - 1L)]
  metadata <- setdiff(metadata, c("Dato", "Transekt"))
  species <- names(all)[seq.int(match("Periode", names(all)) + 1L, ncol(all))]
  all %>% group_by(Flate, Year) %>% summarise(
    across(all_of(setdiff(metadata, c("Flate", "Year", "Blomsterdekke"))),
           ~ {v <- stats::na.omit(.x); if (length(v)) v[[1]] else NA}),
    Blomsterdekke = mean(Blomsterdekke, na.rm = TRUE),
    across(all_of(species), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
}
allBombus_np <- fetch_taxon("Bumblebees")
allLepi_np <- fetch_taxon("Butterflies")
if (!exists("con")) stop("No database connection `con` available after humlesommerfConnect()")
clim_data_raw <- dplyr::tbl(con, DBI::Id(schema = "climate_data", table = "se_norge"))


# ---- Abundance indices -----------------------------------------------------
lepi <- prepare_abundance(allLepi_np, taxon_group = "Butterfly", meta_cols = abund_meta_cols)
bombus <- prepare_abundance(allBombus_np, taxon_group = "Bumblebee", meta_cols = abund_meta_cols)

abund <- bind_rows(lepi, bombus) %>%
  filter(!is.na(year), !is.na(locality_id_abund))

if (is.null(year_min)) year_min <- min(abund$year, na.rm = TRUE)
if (is.null(year_max)) year_max <- max(abund$year, na.rm = TRUE)

abund <- abund %>%
  filter(year >= year_min, year <= year_max)

if (exists("locality_lookup")) {
  locality_lookup_clean <- locality_lookup %>%
    as_tibble() %>%
    clean_names() %>%
    rename(
      locality_id_abund = all_of(janitor::make_clean_names(abund_locality_col)),
      locality_id_clim = all_of(janitor::make_clean_names(clim_locality_col))
    ) %>%
    mutate(
      locality_id_abund = as.character(locality_id_abund),
      locality_id_clim = as.character(locality_id_clim)
    )

  abund <- abund %>%
    left_join(locality_lookup_clean, by = "locality_id_abund")
} else {
  abund <- abund %>%
    mutate(locality_id_clim = locality_id_abund)
}

abund_checks <- abund %>%
  summarise(
    n_rows = n(),
    n_localities = n_distinct(locality_id_abund),
    min_year = min(year),
    max_year = max(year),
    n_zero_total_abundance = sum(total_abundance == 0, na.rm = TRUE)
  )

abund <- abund %>%
  group_by(taxon_group, locality_id_abund) %>%
  mutate(
    log_total_abundance_z = z_score(log1p(total_abundance))
  ) %>%
  ungroup()



# ---- Climate preprocessing -------------------------------------------------
if (inherits(clim_data_raw, "tbl_lazy")) {
  clim_data_raw <- clim_data_raw %>% dplyr::collect()
}

clim <- clim_data_raw %>%
  as_tibble() %>%
  clean_names()

required_climate_cols <- janitor::make_clean_names(c(
  clim_locality_col, clim_date_col, clim_temp_col, clim_prec_col
))

missing_climate_cols <- setdiff(required_climate_cols, names(clim))
if (length(missing_climate_cols) > 0) {
  stop(
    "The following required climate columns are missing after clean_names(): ",
    paste(missing_climate_cols, collapse = ", "),
    ". Available columns are: ",
    paste(names(clim), collapse = ", ")
  )
}

clim <- clim %>%
  rename(
    locality_id_clim = all_of(janitor::make_clean_names(clim_locality_col)),
    date = all_of(janitor::make_clean_names(clim_date_col)),
    temp_mean = all_of(janitor::make_clean_names(clim_temp_col)),
    precip = all_of(janitor::make_clean_names(clim_prec_col))
  ) %>%
  mutate(
    locality_id_clim = as.character(locality_id_clim),
    date = as.Date(date),
    year = lubridate::year(date),
    month = lubridate::month(date),
    doy = lubridate::yday(date),
    temp_mean = as.numeric(temp_mean),
    precip = as.numeric(precip),
    precip = if_else(precip < 0, NA_real_, precip)
  )

clim <- clim %>%
  filter(year >= year_min - 1, year <= year_max)

clim_checks <- clim %>%
  summarise(
    n_rows = n(),
    n_localities = n_distinct(locality_id_clim),
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    missing_temp = sum(is.na(temp_mean)),
    missing_precip = sum(is.na(precip)),
    negative_precip_after_cleaning = sum(precip < 0, na.rm = TRUE)
  )

locality_match_check <- abund %>%
  distinct(locality_id_abund, locality_id_clim) %>%
  mutate(in_climate = locality_id_clim %in% unique(clim$locality_id_clim))


if (any(!locality_match_check$in_climate)) {
  warning(
    "Some abundance localities are missing in climate data. ",
    "Check `locality_match_check`. Abundance uses Flate and climate uses locality; supply `locality_lookup` if necessary."
  )
}


# ---- Seasonal climate windows ----------------------------------------------
clim_winter <- make_climate_window(clim, "winter",
  "01-01",
  "02-28",
  offset_year = 0,
  hot_threshold = 5
)

clim_spring <- make_climate_window(
  clim,
  "spring",
  "03-01",
  "05-31",
  offset_year = 0,
   hot_threshold = 15
)

clim_early <- make_climate_window(
  clim,
  "early",
  "04-01",
  "06-15",
  offset_year = 0,
  hot_threshold = 15
)

clim_summer <- make_climate_window(
  clim,
  "summer",
  "06-01",
  "08-31",
  offset_year = 0,
  hot_threshold = 20
)

clim_growing <- make_climate_window(
  clim,
  "growing",
  "04-01",
  "09-30",
  offset_year = 0,
  hot_threshold = 20
)

clim_prev_late_summer <- make_climate_window(
  clim,
  "prev_late_summer",
  "08-01",
  "09-30",
  offset_year = -1,
  hot_threshold = 18
)

clim_prev_autumn <- make_climate_window(
  clim,
  "prev_autumn",
  "09-01",
  "10-31",
  offset_year = -1,
  hot_threshold = 12
)

clim_summary <- list(clim_winter, clim_spring, clim_early, clim_summer, clim_growing,
                     clim_prev_late_summer, clim_prev_autumn) %>%
                reduce(full_join, by = c("locality_id_clim", "year"))

clim_summary_anom <- add_climate_anomalies(clim_summary, ref_years = anomaly_reference_years)


# ---- Join and standardize climate covariates -------------------------------
dat_clim <- abund %>% left_join(clim_summary_anom, by = c("locality_id_clim", "year")) %>%
  mutate(year_z = z_score(year), year_c = year - mean(year, na.rm = TRUE))
candidate_climate <- c(
  spring_temp = "spring__temp_mean_window__anom__z",
  summer_temp = "summer__temp_mean_window__anom__z",
  spring_precip = "spring__precip_sum__anom__z",
  summer_precip = "summer__precip_sum__anom__z",
  prev_late_temp = "prev_late_summer__temp_mean_window__anom__z",
  prev_late_precip = "prev_late_summer__precip_sum__anom__z")
missing_clim <- setdiff(unname(candidate_climate), names(dat_clim))
if (length(missing_clim)) stop("Missing climate predictors: ", paste(missing_clim, collapse = ", "))
dat_clim <- dat_clim %>% rename(!!!setNames(unname(candidate_climate), names(candidate_climate)))
# Response was standardized within locality and taxon above, as in the original climate model.


# ---- Spatial inputs for the species-level GLLVMs ----------------------------

# Landscape extraction is performed upstream; these variables are used only by
# the two group-level GLLVMs, not by the climate-only glmmTMB model.
spatial_dat <- readRDS(required_file(landscape_path))
if (!"locality" %in% names(spatial_dat)) stop("spatial_dat must contain locality")
landscape <- spatial_dat %>% mutate(locality = as.character(locality))
model_dat <- dat_clim %>%
  left_join(landscape, by = c("locality_id_abund" = "locality"))
names(model_dat) <- gsub(" ", "_", names(model_dat), fixed = TRUE)
if ("habitat_type.x" %in% names(model_dat)) model_dat$habitat_type <- model_dat$habitat_type.x
if (anyDuplicated(model_dat[c("taxon_group", "locality_id_abund", "year")]))
  stop("Duplicate locality-year-taxon rows after landscape join")

# The monitored habitat type is a categorical site-level observation predictor.
# It is separate from land-cover proportions used in GLLVM environmental data.
climate_model_data <- dat_clim %>%
  mutate(
    locality_id_abund = factor(locality_id_abund),
    taxon_group = factor(taxon_group, levels = c("Butterfly", "Bumblebee")),
    habitat_type = factor(habitat_type),
    year = as.integer(year),
    spring_temp2 = spring_temp^2,
    spring_temp3 = spring_temp^3,
    summer_temp2 = summer_temp^2,
    prev_late_temp2 = prev_late_temp^2
  ) %>%
  drop_na(log_total_abundance_z, taxon_group, year, year_z,
          locality_id_abund, habitat_type, spring_temp, summer_temp,
          spring_precip, summer_precip, prev_late_temp, prev_late_precip)

# ---- Final climate-only abundance model -------------------------------------

# Final climate-by-taxon specification (third-order spring temperature),
# including habitat_type and locality-specific temporal AR(1) correlation.
model_abundance_climate <- glmmTMB::glmmTMB(
  log_total_abundance_z ~
    taxon_group * year_z +
    habitat_type +
    taxon_group * spring_temp +
    taxon_group * spring_temp2 +
    taxon_group * spring_temp3 +
    taxon_group * summer_temp +
    taxon_group * summer_temp2 +
    taxon_group * spring_precip +
    taxon_group * summer_precip +
    taxon_group * prev_late_temp +
    taxon_group * prev_late_temp2 +
    taxon_group * prev_late_precip +
    ar1(factor(year) + 0 | locality_id_abund),
  data = climate_model_data,
  family = gaussian()
)

# Land-cover covariates required by the species-level GLLVMs.
land_vars <- function(radius) paste0(c(
  "Settlements_and_other_artificial_areas_", "Forest_and_woodlands_",
  "Forest_and_woodlands_edge_m_", "Cropland_"), radius, "r")


# ---- Species-level GLLVM preparation ---------------------------------------

# Match raw species columns and trait metadata, using 10-observation cutoff.
bb_traits <- readxl::read_excel(required_file(traits_path), sheet = 1)
bf_traits <- readxl::read_excel(required_file(traits_path), sheet = 2)
prepare_species <- function(taxon, trait_df, species_key, prefix = NULL) {
  dat <- model_dat %>% filter(taxon_group == taxon)
  tr <- trait_df %>% mutate(abundance_name = tolower(gsub(" ", "_", as.character(.data[[species_key]]))))
  if (!is.null(prefix)) tr$abundance_name <- sub("^b\\._", "bombus_", tr$abundance_name)
  candidates <- intersect(names(dat), tr$abundance_name)
  if (!length(candidates)) stop("No species columns match the trait workbook for ", taxon)
  y <- as.matrix(dat[, candidates, drop = FALSE]); storage.mode(y) <- "numeric"
  y[is.na(y)] <- 0
  keep <- colSums(y > 0) >= min_species_occurrences
  y <- y[, keep, drop = FALSE]
  if (!ncol(y)) stop("No species left after filtering for ", taxon)
  tr <- tr[match(colnames(y), tr$abundance_name), , drop = FALSE]
  list(dat = dat, y = y, traits = tr)
}
bb <- prepare_species("Bumblebee", bb_traits, "Species", "bombus")
bf <- prepare_species("Butterfly", bf_traits, "Species (not Zygaenidae)")
# GLLVM environmental variables: prior-year climate and landscape at 3000 m.
make_environment <- function(dat) {
  needed <- c("spring_temp", "summer_temp", "prev_late_temp",
              land_vars(3000)[c(1, 4, 2, 3)], "Latitude", "Longitude")
  if (length(setdiff(needed, names(dat)))) stop("Missing GLLVM environmental inputs")
  out <- dat %>% transmute(
    spring_temp, summer_temp, prev_late_summer_temp = prev_late_temp,
    artificial = z_score(.data[[land_vars(3000)[1]]]),
    cropland = z_score(.data[[land_vars(3000)[4]]]),
    forest = z_score(.data[[land_vars(3000)[2]]]),
    forest_edge = z_score(.data[[land_vars(3000)[3]]]),
    latitude = z_score(Latitude), longitude = z_score(Longitude))
  as.data.frame(out)
}
prepare_gllvm_input <- function(dataset, taxon) {
  dat <- dataset$dat
  tr <- dataset$traits
  y <- dataset$y
  required_traits <- if (taxon == "Bumblebee") c("Species temperature index", "Habitat niche width") else
    c("Temperature (STI) - Mean annual temperature in Celsius", "Habitat niche width")
  if (length(setdiff(required_traits, names(tr))))
    stop("Missing expected columns in trait workbook for ", taxon)
  if (taxon == "Bumblebee") {
    tr$STI <- suppressWarnings(as.numeric(tr[["Species temperature index"]]))
    tr$habitat_width <- suppressWarnings(as.numeric(tr[["Habitat niche width"]]))
    # Use a normalized numeric day-of-year trait; otherwise fail explicitly.
    if ("queen_emergence_doy" %in% names(tr))
      tr$phenology <- suppressWarnings(as.numeric(tr$queen_emergence_doy))
    else tr$phenology <- NA_real_
  } else {
    tr$STI <- suppressWarnings(as.numeric(tr[["Temperature (STI) - Mean annual temperature in Celsius"]]))
    tr$habitat_width <- suppressWarnings(as.numeric(tr[["Habitat niche width"]]))
    tr$phenology <- suppressWarnings(as.numeric(tr[["First flight month (Norway)"]]))
  }
  valid <- stats::complete.cases(tr[, c("STI", "habitat_width")])
  y <- y[, valid, drop = FALSE]; tr <- tr[valid, , drop = FALSE]
  env <- make_environment(dat)
  nonempty <- rowSums(y) > 0 & stats::complete.cases(env)
  y <- y[nonempty, , drop = FALSE]
  env <- env[nonempty, , drop = FALSE]
  sites <- factor(dat$locality_id_abund[nonempty])
  trait_mat <- data.frame(STI = z_score(tr$STI), habitat_width = z_score(tr$habitat_width),
                          row.names = colnames(y))
  wide <- data.frame(y, env, site = sites)
  long <- reshape(wide, direction = "long", varying = colnames(y), v.names = "y", timevar = "sp")
  traits_long <- data.frame(sp = seq_len(ncol(y)), species = colnames(y), trait_mat)
  long <- merge(long, traits_long, by = "sp", sort = FALSE)
  long$id <- as.integer(factor(long$id))
  list(long = long, study = data.frame(site = sites), traits = tr, y = y)
}
bb_gllvm <- prepare_gllvm_input(bb, "Bumblebee")
bf_gllvm <- prepare_gllvm_input(bf, "Butterfly")


# ---- Final fourth-corner GLLVMs --------------------------------------------
# Final hypotheses
model_bb_gllvm <- gllvm::gllvm(
  data = bb_gllvm$long,
  formula = y ~ spring_temp + summer_temp + prev_late_summer_temp +
    artificial + cropland + forest + forest_edge + latitude + longitude +
    artificial:habitat_width + cropland:habitat_width,
  family = "negative.binomial", num.lv = 1,
  studyDesign = bb_gllvm$study, row.eff = ~(1 | site),
  randomX = ~ artificial + cropland + forest, seed = 123,
  control.start = list(n.init = 3, jitter.var = 0.01))
model_bf_gllvm <- gllvm::gllvm(
  data = bf_gllvm$long,
  formula = y ~ spring_temp + summer_temp + prev_late_summer_temp +
    artificial + cropland + forest + forest_edge + latitude + longitude +
    spring_temp:STI + summer_temp:STI + prev_late_summer_temp:STI +
    artificial:habitat_width,
  family = "negative.binomial", num.lv = 1,
  studyDesign = bf_gllvm$study, row.eff = ~(1 | site),
  randomX = ~ spring_temp + summer_temp + prev_late_summer_temp +
    artificial + cropland + forest,
  seed = 987, control.start = list(n.init = 10, jitter.var = 0.05))


