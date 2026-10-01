# User settings
settings <- list(
  mode = Sys.getenv("LAAIDS_MODE", "all"),
  data_dir = Sys.getenv("LAAIDS_DATA_DIR", "data"),
  output_dir = Sys.getenv("LAAIDS_OUTPUT_DIR", "outputs"),
  mapping_file = "food_group_mapping.csv",
  bootstrap_reps = as.numeric(Sys.getenv("LAAIDS_BOOT_REPS", "999")),
  mc_reps = as.numeric(Sys.getenv("LAAIDS_MC_REPS", "500")),
  mc_households = as.numeric(Sys.getenv("LAAIDS_MC_HOUSEHOLDS", "600")),
  save_png = TRUE
)

analysis_root <- local({
  sourced_files <- lapply(sys.frames(), function(frame) frame$ofile)
  sourced_files <- Filter(function(x) is.character(x) && length(x) == 1L, sourced_files)
  command_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  file <- if (length(sourced_files)) {
    tail(sourced_files, 1L)[[1L]]
  } else if (length(command_file)) {
    sub("^--file=", "", command_file[1L])
  } else {
    "laaids_public_analysis.R"
  }
  dirname(normalizePath(file, winslash = "/", mustWork = TRUE))
})

# Configuration
resolve_path <- function(path, root) {
  if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", path)) path <- file.path(root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

make_config <- function(settings, root) {
  if (length(settings$mode) != 1L || !settings$mode %in% c("empirical", "monte_carlo", "all")) {
    stop("mode must be 'empirical', 'monte_carlo', or 'all'.", call. = FALSE)
  }
  for (key in c("bootstrap_reps", "mc_reps", "mc_households")) {
    value <- settings[[key]]
    if (length(value) != 1L || is.na(value) || !is.finite(value) || value < 2 || value > .Machine$integer.max || value != floor(value)) {
      stop(key, " must be an integer of at least 2.", call. = FALSE)
    }
    settings[[key]] <- as.integer(value)
  }
  if (settings$mc_households < 10L) stop("mc_households must be at least 10.", call. = FALSE)
  output <- resolve_path(settings$output_dir, root)
  list(
    version = "1.0.0",
    mode = settings$mode,
    data_dir = resolve_path(settings$data_dir, root),
    output_dir = output,
    group_mapping_path = resolve_path(settings$mapping_file, root),
    input_cache_dir = file.path(output, "processed", "input"),
    empirical_dir = file.path(output, "empirical"),
    robustness_dir = file.path(output, "robustness"),
    figure_dir = file.path(output, "figures"),
    audit_dir = file.path(output, "diagnostics", "empirical"),
    mc_output_dir = file.path(output, "monte_carlo"),
    raw_winsor_p = 0.01,
    item_winsor_p = 0.01,
    group_winsor_p = 0.01,
    baseline_min_conversion_coverage = 0.90,
    run_empirical_bootstrap = TRUE,
    empirical_bootstrap_reps = settings$bootstrap_reps,
    empirical_bootstrap_seed = 20260902L,
    empirical_bootstrap_cluster_var = "ea_id",
    empirical_bootstrap_strata_var = "district",
    empirical_bootstrap_require_strata = TRUE,
    empirical_bootstrap_expected_strata = 32L,
    empirical_bootstrap_min_success_rate = 0.90,
    empirical_bootstrap_progress_every = 100L,
    mc_reps = settings$mc_reps,
    mc_households = settings$mc_households,
    mc_numerical_tol = 1e-10,
    mc_scenario_filter = character(0),
    mc_write_raw_csv = FALSE,
    write_local_microdata_csv = FALSE,
    save_png = settings$save_png
  )
}

check_inputs <- function(config) {
  empirical <- config$mode %in% c("empirical", "all")
  packages <- if (empirical) {
    c("haven", "dplyr", "tidyr", "stringr", "lubridate", "readr", "purrr",
      "tibble", "rlang", "broom", "ggplot2", "forcats", "scales")
  } else {
    c("readr", "tibble", "ggplot2")
  }
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) {
    stop("Install required packages before running: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
  if (empirical) {
    files <- c(
      "hh_mod_a_filt.dta", "HH_MOD_META.dta", "HH_MOD_B.dta", "HH_MOD_C.dta",
      "ihs5_consumption_aggregate.dta", "householdgeovariables_ihs5.dta", "HH_MOD_G1.dta",
      "ihs_foodconversion_factor_2020.dta", "mrk_mod_a.dta", "mrk_mod_d.dta", "mrk_mod_otherspecify.dta"
    )
    missing_files <- files[!file.exists(file.path(config$data_dir, files))]
    if (length(missing_files)) {
      stop("Missing input files in ", config$data_dir, ":\n", paste(missing_files, collapse = "\n"), call. = FALSE)
    }
    if (!file.exists(config$group_mapping_path)) stop("Missing food-group mapping: ", config$group_mapping_path, call. = FALSE)
  }
  if (file.exists(config$output_dir) && !dir.exists(config$output_dir)) {
    stop("output_dir is an existing file.", call. = FALSE)
  }
  if (dir.exists(config$output_dir) && length(list.files(config$output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is not empty. Choose a new output_dir to keep runs separate: ", config$output_dir, call. = FALSE)
  }
  suppressPackageStartupMessages(invisible(lapply(packages, library, character.only = TRUE)))
  invisible(TRUE)
}

write_run_metadata <- function(config, status, elapsed_seconds) {
  metadata <- c(config, list(
    run_status = status,
    elapsed_seconds = elapsed_seconds,
    manuscript_replication_counts =
      (config$mode == "monte_carlo" || config$empirical_bootstrap_reps == 999L) &&
      (config$mode == "empirical" || (config$mc_reps == 500L && config$mc_households == 600L))
  ))
  table <- data.frame(
    setting = names(metadata),
    value = vapply(metadata, function(x) paste(x, collapse = ";"), character(1))
  )
  utils::write.csv(table, file.path(config$output_dir, "run_settings.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(config$output_dir, "sessionInfo.txt"))
}

# Input construction
build_input_data <- function(config) {
  data_dir <- config$data_dir
  out_dir <- config$input_cache_dir
  winsor_p <- config$raw_winsor_p
  use_manual_factors <- TRUE
  use_approximate_liquid_factors <- TRUE
  drop_unmatched_conversions <- TRUE
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  required_files <- c(
    "hh_mod_a_filt.dta", "HH_MOD_META.dta", "HH_MOD_B.dta",
    "HH_MOD_C.dta", "ihs5_consumption_aggregate.dta",
    "householdgeovariables_ihs5.dta", "HH_MOD_G1.dta",
    "ihs_foodconversion_factor_2020.dta", "mrk_mod_a.dta",
    "mrk_mod_d.dta", "mrk_mod_otherspecify.dta"
  )
  missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]
  if (length(missing_files)) {
    stop("Required input files are missing from ", data_dir, ": ",
         paste(missing_files, collapse = ", "), call. = FALSE)
  }

  read_dta2 <- function(filename) {
    path <- file.path(data_dir, filename)
    if (!file.exists(path)) stop("File not found: ", path)
    haven::read_dta(path)
  }

  val_chr <- function(x) {
    x <- haven::zap_labels(x)
    if (inherits(x, "Date")) return(as.character(x))
    x <- as.character(x)
    str_trim(x)
  }

  val_num <- function(x) {
    suppressWarnings(as.numeric(haven::zap_labels(x)))
  }

  clean_code <- function(x) {
    x <- val_chr(x)
    x <- str_to_upper(x)
    x <- str_replace_all(x, "\\s+", "")
    x <- na_if(x, "")
    x <- na_if(x, "NA")
    x
  }

  clean_text <- function(x) {
    x <- val_chr(x)
    x <- str_to_upper(x)
    x <- str_replace_all(x, "[^A-Z0-9]+", " ")
    x <- str_squish(x)
    x <- na_if(x, "")
    x <- na_if(x, "NA")
    x
  }

  safe_date <- function(x) {
    if (inherits(x, "Date")) return(x)
    if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as.Date(x))

    z <- suppressWarnings(lubridate::ymd(x))
    if (all(is.na(z))) z <- suppressWarnings(lubridate::dmy(x))
    z
  }

  wtd_mean <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (!any(ok)) return(NA_real_)
    sum(x[ok] * w[ok]) / sum(w[ok])
  }

  wtd_median <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (!any(ok)) return(NA_real_)
    x <- x[ok]; w <- w[ok]
    o <- order(x)
    x <- x[o]; w <- w[o]
    cw <- cumsum(w) / sum(w)
    x[which(cw >= 0.5)[1]]
  }

  winsorize_by <- function(df, value, by, p = 0.01) {

    value_name <- rlang::as_name(rlang::ensym(value))
    by_names <- purrr::map_chr(rlang::ensyms(by), rlang::as_name)

    df %>%
      group_by(across(all_of(by_names))) %>%
      mutate(
        .lo = suppressWarnings(
          quantile(.data[[value_name]], probs = p, na.rm = TRUE, names = FALSE)
        ),
        .hi = suppressWarnings(
          quantile(.data[[value_name]], probs = 1 - p, na.rm = TRUE, names = FALSE)
        ),
        !!value_name := if_else(
          is.finite(.data[[value_name]]),
          pmin(pmax(.data[[value_name]], .lo), .hi),
          .data[[value_name]]
        )
      ) %>%
      ungroup() %>%
      select(-.lo, -.hi)
  }

  subunit_suffix <- function(size_code, flat_heap_code = NULL) {
    s <- clean_code(size_code)
    fh <- clean_code(flat_heap_code)

    out_size <- case_when(
      s %in% c("1", "A", "SMALL") ~ "A",
      s %in% c("2", "B", "MEDIUM") ~ "B",
      s %in% c("3", "C", "LARGE") ~ "C",
      TRUE ~ NA_character_
    )

    out_fh <- case_when(
      fh %in% c("1", "A", "FLAT") ~ "A",
      fh %in% c("2", "B", "HEAPED", "HEAP") ~ "B",
      TRUE ~ NA_character_
    )
    coalesce(out_size, out_fh)
  }

  hh_a <- read_dta2("hh_mod_a_filt.dta") %>%
    transmute(
      case_id,
      HHID,
      ea_id = val_num(ea_id),
      region = val_num(region),
      district = val_num(district),
      reside = val_num(reside),
      urban = as.integer(reside == 1),
      interview_date = safe_date(interviewDate),
      hh_wgt = val_num(hh_wgt),
      hhsize = val_num(hhsize)
    )

  hh_meta <- read_dta2("HH_MOD_META.dta") %>%
    transmute(
      case_id,
      moduleG_start_date = safe_date(moduleG_start_date),
      food_ym = floor_date(moduleG_start_date, "month")
    )

  hh_b <- read_dta2("HH_MOD_B.dta") %>%
    mutate(
      rel_head = val_num(hh_b04),
      sex = val_num(hh_b03),
      age = val_num(hh_b05a)
    )

  head_b <- hh_b %>%
    filter(rel_head == 1) %>%
    transmute(
      case_id,
      head_pid = PID,
      head_female = as.integer(sex == 2),
      head_age = age
    )

  hh_comp <- hh_b %>%
    group_by(case_id) %>%
    summarise(
      n_members = n(),
      n_child_0_5 = sum(age >= 0 & age <= 5, na.rm = TRUE),
      n_child_6_14 = sum(age >= 6 & age <= 14, na.rm = TRUE),
      n_adult_15_64 = sum(age >= 15 & age <= 64, na.rm = TRUE),
      n_elderly_65p = sum(age >= 65, na.rm = TRUE),
      dependency_ratio = (n_child_0_5 + n_child_6_14 + n_elderly_65p) / pmax(n_adult_15_64, 1),
      .groups = "drop"
    )

  hh_c <- read_dta2("HH_MOD_C.dta") %>%
    transmute(
      case_id,
      PID,
      ever_school = val_num(hh_c06),
      highest_grade = val_num(hh_c08),
      highest_qual = val_num(hh_c09),
      can_read = val_num(hh_c05_1),
      can_write = val_num(hh_c05_3)
    )

  head_c <- head_b %>%
    select(case_id, PID = head_pid) %>%
    left_join(hh_c, by = c("case_id", "PID")) %>%
    transmute(
      case_id,
      head_ever_school = ever_school,
      head_highest_grade = highest_grade,
      head_highest_qual = highest_qual,
      head_literate = as.integer(can_read == 1 & can_write == 1)
    )

  consagg <- read_dta2("ihs5_consumption_aggregate.dta") %>%
    transmute(
      case_id,
      adulteq = val_num(adulteq),
      expagg = val_num(expagg),
      rexpagg = val_num(rexpagg),
      expaggpc = val_num(expaggpc),
      rexpaggpc = val_num(rexpaggpc),
      poor = val_num(poor)
    )

  geo <- read_dta2("householdgeovariables_ihs5.dta") %>%
    transmute(
      case_id,
      dist_road = val_num(dist_road),
      dist_agmrkt = val_num(dist_agmrkt),
      dist_admarc = val_num(dist_admarc),
      dist_popcenter = val_num(dist_popcenter),
      dist_boma = val_num(dist_boma)
    )

  hh_base <- hh_a %>%
    left_join(hh_meta, by = "case_id") %>%
    left_join(head_b, by = "case_id") %>%
    left_join(head_c, by = "case_id") %>%
    left_join(hh_comp, by = "case_id") %>%
    left_join(consagg, by = "case_id") %>%
    left_join(geo, by = "case_id") %>%
    mutate(
      food_ym = coalesce(food_ym, floor_date(interview_date, "month")),
      food_year = year(food_ym),
      food_month = month(food_ym)
    )

  readr::write_csv(hh_base, file.path(out_dir, "hh_base_controls.csv"))

  food_raw <- read_dta2("HH_MOD_G1.dta") %>%
    mutate(
      consumed = val_num(hh_g01),
      item_code = val_num(hh_g02),
      item_name = as.character(haven::as_factor(hh_g02, levels = "labels")),
      q_purch_raw = val_num(hh_g04a),
      unit_code_original = clean_code(hh_g04b),
      unit_label = val_chr(hh_g04b_label),
      unit_other = clean_text(hh_g04b_oth),
      sub_size = clean_code(hh_g04c),
      sub_flatheap = clean_code(hh_g04c_1),
      exp_purch = val_num(hh_g05),
      unit_suffix = subunit_suffix(hh_g04c, hh_g04c_1),
      unit_code_alt = if_else(
        !is.na(unit_suffix) & !str_detect(unit_code_original, "[A-Z]$"),
        paste0(unit_code_original, unit_suffix),
        unit_code_original
      )
    ) %>%
    filter(consumed == 1, q_purch_raw > 0, exp_purch > 0) %>%
    left_join(hh_base %>% select(case_id, region, district, ea_id, food_ym, hh_wgt), by = "case_id") %>%
    mutate(
      uv_raw = exp_purch / q_purch_raw,
      log_uv_raw = log(uv_raw)
    )

  item_dict <- food_raw %>%
    distinct(item_code, item_name) %>%
    arrange(item_code)
  readr::write_csv(item_dict, file.path(out_dir, "food_item_dictionary.csv"))

  conv <- read_dta2("ihs_foodconversion_factor_2020.dta") %>%
    transmute(
      region = val_num(region),
      item_code = val_num(item_code),
      item_name_conv = val_chr(item_name),
      unit_name_conv = val_chr(unit_name),
      unit_code = clean_code(unit_code),
      other_unit = clean_text(Otherunit),
      factor = val_num(factor),
      source = val_num(source)
    ) %>%
    filter(!is.na(region), !is.na(item_code), !is.na(unit_code), is.finite(factor), factor > 0)

  conv_nonother <- conv %>%
    filter(is.na(other_unit)) %>%
    group_by(region, item_code, unit_code) %>%
    summarise(
      factor = median(factor, na.rm = TRUE),
      source = suppressWarnings(min(source, na.rm = TRUE)),
      .groups = "drop"
    )

  conv_other <- conv %>%
    filter(!is.na(other_unit)) %>%
    group_by(region, item_code, unit_code, other_unit) %>%
    summarise(
      factor = median(factor, na.rm = TRUE),
      source = suppressWarnings(min(source, na.rm = TRUE)),
      .groups = "drop"
    )

  food_conv1 <- food_raw %>%
    left_join(
      conv_nonother %>% rename(factor_1 = factor, source_1 = source),
      by = c("region", "item_code", "unit_code_original" = "unit_code")
    )

  food_conv2 <- food_conv1 %>%
    left_join(
      conv_nonother %>% rename(factor_2 = factor, source_2 = source),
      by = c("region", "item_code", "unit_code_alt" = "unit_code")
    )

  food_conv3 <- food_conv2 %>%
    left_join(
      conv_other %>% rename(factor_3 = factor, source_3 = source),
      by = c("region", "item_code", "unit_code_original" = "unit_code", "unit_other" = "other_unit")
    )

  food_conv <- food_conv3 %>%
    mutate(
      item_num = suppressWarnings(as.integer(item_code)),
      factor_manual_strict = case_when(
        unit_code_original %in% c("1")  ~ 1,
        unit_code_original %in% c("18") ~ 0.001,
        TRUE ~ NA_real_
      ),
      factor_manual_alimi = case_when(

        unit_code_original %in% c("1")  ~ 1,
        unit_code_original %in% c("18") ~ 0.001,

        use_manual_factors & unit_code_original %in% c("31") ~ 0.3,
        use_manual_factors & unit_code_original %in% c("32") ~ 0.6,
        use_manual_factors & unit_code_original %in% c("33") ~ 0.7,
        use_manual_factors & unit_code_original %in% c("34") ~ 0.15,
        use_manual_factors & unit_code_original %in% c("35") ~ 0.4,
        use_manual_factors & unit_code_original %in% c("36") ~ 0.5,
        use_manual_factors & unit_code_original %in% c("37") ~ 1,
        use_manual_factors & unit_code_original %in% c("41") ~ 0.025,
        use_manual_factors & unit_code_original %in% c("42") ~ 0.05,
        use_manual_factors & unit_code_original %in% c("43") ~ 0.1,
        use_manual_factors & unit_code_original %in% c("65") ~ 0.25,
        use_manual_factors & unit_code_original %in% c("70") ~ 0.025,
        use_manual_factors & unit_code_original %in% c("71") ~ 0.1,
        use_manual_factors & unit_code_original %in% c("72") ~ 0.25,
        use_manual_factors & unit_code_original %in% c("73") ~ 0.5,

        use_manual_factors & unit_code_original == "15" & (item_num %in% c(705, 706, 813) | dplyr::between(item_num, 904, 915)) ~ 1,
        use_manual_factors & unit_code_original == "19" & (item_num %in% c(705, 706, 813) | dplyr::between(item_num, 904, 915)) ~ 0.001,
        use_manual_factors & unit_code_original == "15" & item_num %in% c(815, 817) ~ 1.4,
        use_manual_factors & unit_code_original == "19" & item_num %in% c(815, 817) ~ 0.001 * 1.4,
        use_manual_factors & unit_code_original == "15" & item_num %in% c(701) ~ 1.1,
        use_manual_factors & unit_code_original == "19" & item_num %in% c(701) ~ 0.001 * 1.1,
        use_manual_factors & unit_code_original == "15" & item_num %in% c(803) ~ 0.92,
        use_manual_factors & unit_code_original == "19" & item_num %in% c(803) ~ 0.001 * 0.92,
        use_manual_factors & unit_code_original == "15" & item_num %in% c(814) ~ 1.08,
        use_manual_factors & unit_code_original == "19" & item_num %in% c(814) ~ 0.001 * 1.08,

        use_approximate_liquid_factors & unit_code_original %in% c("15") ~ 1,
        use_approximate_liquid_factors & unit_code_original %in% c("19") ~ 0.001,
        TRUE ~ NA_real_
      ),
      factor_manual = coalesce(factor_manual_alimi, factor_manual_strict),
      factor = coalesce(factor_1, factor_2, factor_3, factor_manual),
      factor_source = case_when(
        !is.na(factor_1) ~ "conversion_direct",
        is.na(factor_1) & !is.na(factor_2) ~ "conversion_subunit_alt",
        is.na(factor_1) & is.na(factor_2) & !is.na(factor_3) ~ "conversion_other_text",
        is.na(factor_1) & is.na(factor_2) & is.na(factor_3) & !is.na(factor_manual) & !is.na(factor_manual_alimi) ~ "manual_alimi_rule",
        is.na(factor_1) & is.na(factor_2) & is.na(factor_3) & !is.na(factor_manual) ~ "manual_standard_unit",
        TRUE ~ "unmatched"
      ),
      q_std = q_purch_raw * factor,
      uv_std = exp_purch / q_std,
      log_uv_std = log(uv_std),
      log_unit_scale = log_uv_raw - log_uv_std,
      conversion_matched = is.finite(factor) & factor > 0
    )

  conversion_coverage <- food_conv %>%
    summarise(
      purchase_rows = n(),
      matched_rows = sum(conversion_matched, na.rm = TRUE),
      matched_row_share = matched_rows / purchase_rows,
      purchase_exp_total = sum(exp_purch, na.rm = TRUE),
      purchase_exp_matched = sum(exp_purch[conversion_matched], na.rm = TRUE),
      matched_expenditure_share = purchase_exp_matched / purchase_exp_total,
      unmatched_rows = sum(!conversion_matched, na.rm = TRUE),
      unmatched_expenditure = sum(exp_purch[!conversion_matched], na.rm = TRUE)
    )
  readr::write_csv(conversion_coverage, file.path(out_dir, "conversion_coverage_summary.csv"))

  unmatched_units <- food_conv %>%
    filter(!conversion_matched) %>%
    group_by(item_code, item_name, unit_code_original, unit_label, unit_other, unit_code_alt) %>%
    summarise(
      n = n(),
      exp = sum(exp_purch, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(exp))
  readr::write_csv(unmatched_units, file.path(out_dir, "unmatched_conversion_units.csv"))

  conversion_sources <- food_conv %>%
    group_by(factor_source) %>%
    summarise(
      n = n(),
      exp = sum(exp_purch, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(row_share = n / sum(n), exp_share = exp / sum(exp)) %>%
    arrange(desc(exp))
  readr::write_csv(conversion_sources, file.path(out_dir, "conversion_factor_sources.csv"))

  hh_conversion_coverage <- food_conv %>%
    group_by(case_id) %>%
    summarise(
      total_food_purchase_all = sum(exp_purch, na.rm = TRUE),
      total_food_purchase_matched = sum(exp_purch[conversion_matched], na.rm = TRUE),
      matched_exp_share_h = total_food_purchase_matched / total_food_purchase_all,
      n_items_purchased_all = n_distinct(item_code),
      n_items_conversion_matched = n_distinct(item_code[conversion_matched]),
      .groups = "drop"
    )
  readr::write_csv(hh_conversion_coverage, file.path(out_dir, "household_conversion_coverage.csv"))

  if (drop_unmatched_conversions) {
    food_uv <- food_conv %>% filter(conversion_matched)
  } else {
    food_uv <- food_conv
  }

  food_uv <- food_uv %>%
    filter(is.finite(uv_raw), uv_raw > 0, is.finite(uv_std), uv_std > 0) %>%
    winsorize_by(uv_raw, by = item_code, p = winsor_p) %>%
    winsorize_by(uv_std, by = item_code, p = winsor_p) %>%
    mutate(
      log_uv_raw = log(uv_raw),
      log_uv_std = log(uv_std),
      log_unit_scale = log_uv_raw - log_uv_std
    )

  food_uv <- food_uv %>%
    group_by(case_id) %>%
    mutate(
      total_food_purchase = sum(exp_purch, na.rm = TRUE),
      w_item = exp_purch / total_food_purchase
    ) %>%
    ungroup() %>%
    filter(is.finite(w_item), w_item > 0)

  market_a <- read_dta2("mrk_mod_a.dta") %>%
    transmute(
      MarketID = val_num(MarketID),
      MarketID_Visit = val_chr(MarketID_Visit),
      mrk_region = val_num(mrk_region),
      dist_name = val_chr(dist_name),
      reside_mrk = val_num(reside),
      market_date = safe_date(mrk_a08),
      market_month_var = val_num(month),
      market_ym = floor_date(market_date, "month")
    )

  market_a <- market_a %>%
    mutate(
      market_ym = if_else(
        is.na(market_ym) & is.finite(market_month_var),
        as.Date(sprintf("2019-%02d-01", pmax(pmin(as.integer(market_month_var), 12L), 1L))),
        market_ym
      )
    )

  map_market_item_to_hh <- function(code) {
    code_chr <- str_to_upper(str_trim(as.character(code)))

    manual <- c(
      "119" = "102", "120" = "102",
      "316A" = "301", "316B" = "301", "316C" = "301",
      "806" = "811", "807" = "811", "808" = "811",
      "850" = "803", "851" = "803", "852" = "803"
    )
    out <- ifelse(code_chr %in% names(manual), manual[code_chr], NA_character_)

    numeric_part <- str_extract(code_chr, "^[0-9]+")
    out <- ifelse(is.na(out), numeric_part, out)
    suppressWarnings(as.integer(out))
  }

  make_market_prices_from_d <- function(df, source_name = "D") {
    df %>%
      transmute(
        MarketID = val_num(MarketID),
        MarketID_Visit = val_chr(MarketID_Visit),
        avail = val_num(D01),
        item_code_market = str_to_upper(str_trim(val_chr(d1_b))),
        item_name_market = val_chr(d1_a),
        unit_name_market = val_chr(d1_c),
        unit_code_market = clean_code(d1_d),
        kg1 = val_num(d1_g), price1 = val_num(d1_h),
        kg2 = val_num(d1_i), price2 = val_num(d1_j),
        kg3 = val_num(d1_k), price3 = val_num(d1_l),
        readily_available = val_num(d1_m),
        source_module = source_name
      ) %>%
      mutate(
        item_code = map_market_item_to_hh(item_code_market),
        p1 = if_else(is.finite(kg1) & kg1 > 0 & is.finite(price1) & price1 > 0, price1 / kg1, NA_real_),
        p2 = if_else(is.finite(kg2) & kg2 > 0 & is.finite(price2) & price2 > 0, price2 / kg2, NA_real_),
        p3 = if_else(is.finite(kg3) & kg3 > 0 & is.finite(price3) & price3 > 0, price3 / kg3, NA_real_),
        price_kg_row = pmap_dbl(list(p1, p2, p3), ~ median(c(...), na.rm = TRUE)),
        price_kg_row = if_else(is.nan(price_kg_row), NA_real_, price_kg_row)
      ) %>%
      filter(!is.na(item_code), is.finite(price_kg_row), price_kg_row > 0)
  }

  mrk_d <- read_dta2("mrk_mod_d.dta")
  market_prices_d <- make_market_prices_from_d(mrk_d, "D")

  other_path <- file.path(data_dir, "mrk_mod_otherspecify.dta")
  if (file.exists(other_path)) {
    mrk_other <- read_dta2("mrk_mod_otherspecify.dta")

    if (!"D01" %in% names(mrk_other)) mrk_other$D01 <- NA_real_
    market_prices_other <- make_market_prices_from_d(mrk_other, "D_other")
    market_prices_row <- bind_rows(market_prices_d, market_prices_other)
  } else {
    market_prices_row <- market_prices_d
  }

  market_prices_row <- market_prices_row %>%
    left_join(market_a, by = c("MarketID", "MarketID_Visit")) %>%
    filter(!is.na(market_ym)) %>%
    winsorize_by(price_kg_row, by = item_code, p = winsor_p) %>%
    mutate(log_market_price = log(price_kg_row))

  market_mapping_review <- market_prices_row %>%
    distinct(item_code_market, item_code, item_name_market, unit_name_market, unit_code_market, source_module) %>%
    arrange(item_code, item_code_market, unit_code_market)
  readr::write_csv(market_mapping_review, file.path(out_dir, "market_item_mapping_review.csv"))

  market_price_region_month <- market_prices_row %>%
    group_by(item_code, mrk_region, market_ym) %>%
    summarise(
      market_price_rm = median(price_kg_row, na.rm = TRUE),
      n_market_obs_rm = n(),
      .groups = "drop"
    ) %>%
    mutate(log_market_price_rm = log(market_price_rm))

  market_price_nat_month <- market_prices_row %>%
    group_by(item_code, market_ym) %>%
    summarise(
      market_price_nm = median(price_kg_row, na.rm = TRUE),
      n_market_obs_nm = n(),
      .groups = "drop"
    ) %>%
    mutate(log_market_price_nm = log(market_price_nm))

  market_price_item_nat <- market_prices_row %>%
    group_by(item_code) %>%
    summarise(
      market_price_n = median(price_kg_row, na.rm = TRUE),
      n_market_obs_n = n(),
      .groups = "drop"
    ) %>%
    mutate(log_market_price_n = log(market_price_n))

  saveRDS(market_prices_row, file.path(out_dir, "market_prices_vendor_median_rows.rds"))
  readr::write_csv(market_price_region_month, file.path(out_dir, "market_price_item_region_month.csv"))
  readr::write_csv(market_price_nat_month, file.path(out_dir, "market_price_item_national_month.csv"))
  readr::write_csv(market_price_item_nat, file.path(out_dir, "market_price_item_national_allmonths.csv"))

  market_summary <- tibble(
    n_price_rows = nrow(market_prices_row),
    n_items = n_distinct(market_prices_row$item_code),
    n_market_visits = n_distinct(market_prices_row$MarketID_Visit),
    n_markets = n_distinct(market_prices_row$MarketID)
  )
  readr::write_csv(market_summary, file.path(out_dir, "market_price_summary.csv"))

  analysis_item <- food_uv %>%
    left_join(
      market_price_region_month,
      by = c("item_code", "region" = "mrk_region", "food_ym" = "market_ym")
    ) %>%
    left_join(
      market_price_nat_month,
      by = c("item_code", "food_ym" = "market_ym")
    ) %>%
    left_join(market_price_item_nat, by = "item_code") %>%
    mutate(
      log_market_price = coalesce(log_market_price_rm, log_market_price_nm, log_market_price_n),
      market_price_match_level = case_when(
        !is.na(log_market_price_rm) ~ "region_month",
        is.na(log_market_price_rm) & !is.na(log_market_price_nm) ~ "national_month",
        is.na(log_market_price_rm) & is.na(log_market_price_nm) & !is.na(log_market_price_n) ~ "national_allmonths",
        TRUE ~ "unmatched"
      )
    )

  price_match_summary <- analysis_item %>%
    group_by(market_price_match_level) %>%
    summarise(
      n = n(),
      exp = sum(exp_purch, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(row_share = n / sum(n), exp_share = exp / sum(exp))
  readr::write_csv(price_match_summary, file.path(out_dir, "market_price_match_summary.csv"))

  saveRDS(analysis_item, file.path(out_dir, "analysis_item_level_uv_market_prices.rds"))
  invisible(out_dir)
}

# Measurement units
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
  invisible(path)
}

require_columns <- function(data, columns, object_name = deparse(substitute(data))) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      object_name, " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

add_missing_columns <- function(data, columns, value = NA_real_) {
  for (column in setdiff(columns, names(data))) data[[column]] <- value
  data
}

finite <- function(x) x[is.finite(x)]

safe_mean <- function(x) {
  x <- finite(x)
  if (length(x) == 0) NA_real_ else mean(x)
}

safe_median <- function(x) {
  x <- finite(x)
  if (length(x) == 0) NA_real_ else stats::median(x)
}

safe_quantile <- function(x, probability) {
  x <- finite(x)
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, probability, names = FALSE, na.rm = TRUE))
}

wtd_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

wtd_cor <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  if (sum(ok) < 3) return(NA_real_)
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  w <- w / sum(w)
  mx <- sum(w * x)
  my <- sum(w * y)
  vx <- sum(w * (x - mx)^2)
  vy <- sum(w * (y - my)^2)
  if (!is.finite(vx) || !is.finite(vy) || vx <= 0 || vy <= 0) return(NA_real_)
  sum(w * (x - mx) * (y - my)) / sqrt(vx * vy)
}

wtd_mae <- function(x, y, w) {
  wtd_mean(abs(x - y), w)
}

wtd_rmse <- function(x, y, w) {
  value <- wtd_mean((x - y)^2, w)
  if (!is.finite(value)) NA_real_ else sqrt(value)
}

winsorize_by <- function(data, value, by, p = 0.01, suffix = NULL) {
  value_name <- rlang::as_name(rlang::ensym(value))
  by_names <- purrr::map_chr(rlang::ensyms(by), rlang::as_name)
  target_name <- if (is.null(suffix)) value_name else paste0(value_name, suffix)

  data %>%
    group_by(across(all_of(by_names))) %>%
    mutate(
      .lo = safe_quantile(.data[[value_name]], p),
      .hi = safe_quantile(.data[[value_name]], 1 - p),
      !!target_name := case_when(
        !is.finite(.data[[value_name]]) ~ .data[[value_name]],
        !is.finite(.lo) | !is.finite(.hi) ~ .data[[value_name]],
        TRUE ~ pmin(pmax(.data[[value_name]], .lo), .hi)
      )
    ) %>%
    ungroup() %>%
    select(-.lo, -.hi)
}

winsor_thresholds_by <- function(data, value, by, p = 0.01, stage) {
  value_name <- rlang::as_name(rlang::ensym(value))
  by_names <- purrr::map_chr(rlang::ensyms(by), rlang::as_name)

  data %>%
    group_by(across(all_of(by_names))) %>%
    summarise(
      stage = stage,
      source_variable = value_name,
      winsor_probability_each_tail = p,
      n_rows = n(),
      n_finite = sum(is.finite(.data[[value_name]])),
      lower_threshold = safe_quantile(.data[[value_name]], p),
      upper_threshold = safe_quantile(.data[[value_name]], 1 - p),
      .groups = "drop"
    )
}

assert_unique_key <- function(data, key, object_name) {
  require_columns(data, key, object_name)
  invalid_key <- rep(FALSE, nrow(data))
  for (column in key) {
    values <- data[[column]]
    invalid_key <- invalid_key | is.na(values)
    if (is.character(values)) {
      invalid_key <- invalid_key | trimws(values) == ""
    }
  }
  if (any(invalid_key)) {
    stop(
      object_name, " contains ", sum(invalid_key),
      " row(s) with a missing or blank key in ",
      paste(key, collapse = " + "), ".",
      call. = FALSE
    )
  }

  duplicate_rows <- data %>%
    group_by(across(all_of(key))) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)
  if (nrow(duplicate_rows) > 0) {
    stop(
      object_name, " must contain at most one row per ",
      paste(key, collapse = " + "), ". Found ", nrow(duplicate_rows),
      " duplicated key(s).",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

clean_unit_code <- function(x) {
  x <- as.character(x)
  x <- stringr::str_to_upper(stringr::str_trim(x))
  x <- stringr::str_replace_all(x, "\\s+", "")
  x[x %in% c("", "NA", "N/A")] <- NA_character_

  numeric_like <- !is.na(x) & stringr::str_detect(x, "^[0-9]+(?:\\.0+)?$")
  numeric_code <- suppressWarnings(as.numeric(x[numeric_like]))
  x[numeric_like] <- ifelse(
    is.finite(numeric_code),
    format(as.integer(numeric_code), scientific = FALSE, trim = TRUE),
    x[numeric_like]
  )
  x
}

classify_factor_source <- function(data) {
  data <- add_missing_columns(
    data,
    c("factor_1", "factor_2", "factor_3", "factor_manual_alimi", "factor_manual_strict")
  )
  if (!"factor_source" %in% names(data)) data$factor_source <- NA_character_
  for (column in c("factor_1", "factor_2", "factor_3")) {
    data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
  }

  package_codes <- c("31", "32", "33", "34", "35", "36", "37", "41", "42", "43", "65", "70", "71", "72", "73")
  liquid_codes <- c("15", "19")
  item_code_num <- suppressWarnings(as.integer(data$item_code))
  item_specific_liquid <- (
    item_code_num %in% c(705, 706, 813, 815, 817, 701, 803, 814) |
      dplyr::between(item_code_num, 904L, 915L)
  )

  data %>%
    mutate(
      unit_code_original = clean_unit_code(unit_code_original),
      factor_source_detailed = case_when(
        is.finite(factor_1) & factor_1 > 0 ~ "official_direct",
        (!is.finite(factor_1) | factor_1 <= 0) & is.finite(factor_2) & factor_2 > 0 ~ "official_subunit_code",
        (!is.finite(factor_1) | factor_1 <= 0) &
          (!is.finite(factor_2) | factor_2 <= 0) &
          is.finite(factor_3) & factor_3 > 0 ~ "official_other_text",
        unit_code_original == "1" ~ "exact_kilogram",
        unit_code_original == "18" ~ "exact_gram",
        unit_code_original %in% package_codes ~ "manual_package_or_tin",
        unit_code_original %in% liquid_codes & item_specific_liquid ~ "manual_item_specific_liquid",
        unit_code_original %in% liquid_codes ~ "rough_liquid_fallback",
        factor_source == "manual_standard_unit" ~ "other_standard_unit_rule",
        factor_source == "manual_alimi_rule" ~ "other_manual_rule",
        !is.na(factor_source) ~ as.character(factor_source),
        TRUE ~ "unclassified_matched"
      ),
      factor_source_is_official_or_exact = factor_source_detailed %in% c(
        "official_direct", "official_subunit_code", "official_other_text",
        "exact_kilogram", "exact_gram"
      ),
      factor_source_is_manual_liquid = factor_source_detailed %in% c(
        "manual_item_specific_liquid", "rough_liquid_fallback"
      )
    )
}

prepare_paired_unit_values <- function(data, winsor_p = 0.01) {
  require_columns(
    data,
    c(
      "case_id", "item_code", "exp_purch", "q_purch_raw", "factor",
      "unit_code_original", "log_market_price"
    ),
    "analysis_item"
  )

  data <- classify_factor_source(data)
  if (!"log_uv_raw" %in% names(data)) data$log_uv_raw <- NA_real_
  if (!"log_uv_std" %in% names(data)) data$log_uv_std <- NA_real_

  data <- data %>%
    mutate(
      case_id = as.character(case_id),
      item_code = as.integer(item_code),
      factor = as.numeric(factor),
      exp_purch = as.numeric(exp_purch),
      q_purch_raw = as.numeric(q_purch_raw),
      uv_raw_original = exp_purch / q_purch_raw,
      q_std_original = q_purch_raw * factor,
      uv_std_original = exp_purch / q_std_original,
      log_uv_raw_original = log(uv_raw_original),
      log_uv_std_original = log(uv_std_original),
      log_factor = log(factor)
    ) %>%
    filter(
      is.finite(exp_purch), exp_purch > 0,
      is.finite(q_purch_raw), q_purch_raw > 0,
      is.finite(factor), factor > 0,
      is.finite(log_uv_std_original), is.finite(log_factor)
    )

  item_winsorization_thresholds <- winsor_thresholds_by(
    data,
    log_uv_std_original,
    by = item_code,
    p = winsor_p,
    stage = "item_level_common_unit_log_value"
  )

  data <- winsorize_by(
    data,
    log_uv_std_original,
    by = item_code,
    p = winsor_p,
    suffix = "_paired_clean"
  ) %>%
    mutate(

      log_uv_std = log_uv_std_original_paired_clean,
      log_uv_raw = log_uv_std + log_factor,
      uv_std = exp(log_uv_std),
      uv_raw = exp(log_uv_raw),
      q_std = exp_purch / uv_std,
      paired_cleaning_changed = abs(log_uv_std - log_uv_std_original) > 1e-12,
      paired_identity_error = (log_uv_raw - log_uv_std) - log_factor,
      is_nonkg_reported_unit = !is.na(unit_code_original) & unit_code_original != "1"
    ) %>%
    select(-log_uv_std_original_paired_clean)

  max_error <- max(abs(data$paired_identity_error), na.rm = TRUE)
  if (!is.finite(max_error) || max_error > 1e-10) {
    stop(
      "Paired unit-value identity failed after reconstruction. Maximum absolute error = ",
      format(max_error, scientific = TRUE),
      call. = FALSE
    )
  }

  attr(data, "item_winsorization_thresholds") <- item_winsorization_thresholds
  data
}

compute_household_exposure <- function(item_data, coverage_data, hh_base) {
  require_columns(
    item_data,
    c("case_id", "item_code", "exp_purch", "log_factor", "is_nonkg_reported_unit"),
    "item_data"
  )
  require_columns(
    coverage_data,
    c("case_id", "total_food_purchase_all"),
    "coverage_data"
  )
  require_columns(hh_base, "case_id", "hh_base")

  coverage_data <- coverage_data %>% mutate(case_id = as.character(case_id))
  hh_base <- hh_base %>% mutate(case_id = as.character(case_id))
  assert_unique_key(coverage_data, "case_id", "coverage_data")
  assert_unique_key(hh_base, "case_id", "hh_base")

  if (!"hh_wgt" %in% names(hh_base)) {
    warning(
      "hh_base has no hh_wgt column; equal household weights will be used for weighted summaries.",
      call. = FALSE
    )
    hh_base$hh_wgt <- 1
  }

  household <- item_data %>%
    group_by(case_id) %>%
    summarise(
      total_food_purchase_clean = sum(exp_purch, na.rm = TRUE),
      C_h = sum(exp_purch * log_factor, na.rm = TRUE) / total_food_purchase_clean,
      share_nonkg_exp_matched = sum(exp_purch[is_nonkg_reported_unit], na.rm = TRUE) /
        total_food_purchase_clean,
      n_items_purchased_matched = n_distinct(item_code),
      n_rows_matched = n(),
      .groups = "drop"
    ) %>%
    left_join(coverage_data, by = "case_id") %>%
    mutate(
      total_food_purchase_all = coalesce(
        as.numeric(total_food_purchase_all),
        as.numeric(total_food_purchase_clean)
      ),
      matched_exp_share_clean = total_food_purchase_clean / total_food_purchase_all,
      share_nonkg_exp_all = NA_real_,
      share_nonkg_exp = share_nonkg_exp_matched
    )

  household$share_nonkg_exp_all <- household$share_nonkg_exp_matched *
    household$total_food_purchase_clean / household$total_food_purchase_all

  hh_columns_to_join <- setdiff(
    names(hh_base),
    setdiff(names(household), "case_id")
  )
  household <- household %>%
    left_join(
      hh_base %>% select(all_of(unique(c("case_id", hh_columns_to_join)))),
      by = "case_id"
    )

  C_mean <- mean(household$C_h, na.rm = TRUE)
  C_weighted_mean <- wtd_mean(household$C_h, household$hh_wgt)
  if (!is.finite(C_weighted_mean)) C_weighted_mean <- C_mean

  household %>%
    mutate(
      C_h_centered = C_h - C_mean,
      C_h_centered_weighted = C_h - C_weighted_mean,
      abs_C_h_centered = abs(C_h_centered),
      abs_C_h_centered_weighted = abs(C_h_centered_weighted)
    )
}

rebuild_measurement_units <- function(config) {
  item_path <- file.path(config$input_cache_dir, "analysis_item_level_uv_market_prices.rds")
  coverage_path <- file.path(config$input_cache_dir, "household_conversion_coverage.csv")
  hh_base_path <- file.path(config$input_cache_dir, "hh_base_controls.csv")
  item_dict_path <- file.path(config$input_cache_dir, "food_item_dictionary.csv")

  invisible(lapply(
    c(item_path, coverage_path, hh_base_path, item_dict_path),
    stop_if_missing
  ))

  input_item <- readRDS(item_path)
  household_coverage <- readr::read_csv(coverage_path, show_col_types = FALSE)
  household_base <- readr::read_csv(hh_base_path, show_col_types = FALSE)
  item_dictionary <- readr::read_csv(item_dict_path, show_col_types = FALSE)

  item_data <- prepare_paired_unit_values(
    input_item,
    winsor_p = config$item_winsor_p
  )
  item_winsorization_thresholds <- attr(
    item_data,
    "item_winsorization_thresholds"
  )
  if (is.null(item_winsorization_thresholds)) {
    stop("Item-level winsorization thresholds were not retained.", call. = FALSE)
  }
  attr(item_data, "item_winsorization_thresholds") <- NULL

  household_data <- compute_household_exposure(
    item_data = item_data,
    coverage_data = household_coverage,
    hh_base = household_base
  )

  factor_source_summary <- item_data %>%
    group_by(factor_source_detailed) %>%
    summarise(
      n_rows = n(),
      expenditure = sum(exp_purch, na.rm = TRUE),
      n_households = n_distinct(case_id),
      .groups = "drop"
    ) %>%
    mutate(
      row_share = n_rows / sum(n_rows),
      expenditure_share = expenditure / sum(expenditure)
    ) %>%
    arrange(desc(expenditure))

  conversion_rule_examples <- item_data %>%
    mutate(
      item_name = if ("item_name" %in% names(item_data)) as.character(item_name) else NA_character_,
      unit_label = if ("unit_label" %in% names(item_data)) as.character(unit_label) else NA_character_
    ) %>%
    distinct(
      factor_source_detailed, item_code, item_name,
      unit_code_original, unit_label, factor
    ) %>%
    arrange(factor_source_detailed, item_code, unit_code_original, factor) %>%
    group_by(factor_source_detailed) %>%
    slice_head(n = 5) %>%
    ungroup()

  market_match_source <- item_data
  if ("market_price_match_level" %in% names(market_match_source)) {
    market_match_source$market_price_match_level <- as.character(market_match_source$market_price_match_level)
  } else {
    market_match_source$market_price_match_level <- ifelse(
      is.finite(market_match_source$log_market_price),
      "matched_unspecified", "unmatched"
    )
  }
  market_match_all <- market_match_source %>%
    group_by(market_price_match_level) %>%
    summarise(
      n_rows = n(),
      expenditure = sum(exp_purch, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      row_share = n_rows / sum(n_rows),
      expenditure_share = expenditure / sum(expenditure)
    )

  benchmark_data <- item_data %>%
    filter(
      is.finite(log_market_price),
      is.finite(log_uv_raw),
      is.finite(log_uv_std),
      is.finite(exp_purch), exp_purch > 0
    )

  unit_value_market_benchmark <- purrr::map_dfr(
    c(raw = "log_uv_raw", standardized = "log_uv_std"),
    function(variable) {
      x <- benchmark_data[[variable]]
      y <- benchmark_data$log_market_price
      w <- benchmark_data$exp_purch
      tibble(
        benchmark_source = "separately_collected_market_survey_price",
        unit_value_measure = if (identical(variable, "log_uv_raw")) "raw" else "standardized",
        n_item_observations = sum(is.finite(x) & is.finite(y)),
        expenditure_weighted_correlation = wtd_cor(x, y, w),
        expenditure_weighted_mae = wtd_mae(x, y, w),
        expenditure_weighted_rmse = wtd_rmse(x, y, w),
        unweighted_correlation = suppressWarnings(cor(x, y, use = "complete.obs")),
        unweighted_mae = mean(abs(x - y), na.rm = TRUE),
        unweighted_rmse = sqrt(mean((x - y)^2, na.rm = TRUE))
      )
    }
  )

  benchmark_by_match_source <- benchmark_data
  if ("market_price_match_level" %in% names(benchmark_by_match_source)) {
    benchmark_by_match_source$market_price_match_level <- as.character(benchmark_by_match_source$market_price_match_level)
  } else {
    benchmark_by_match_source$market_price_match_level <- "matched_unspecified"
  }
  benchmark_by_match_level <- benchmark_by_match_source %>%
    group_by(market_price_match_level) %>%
    group_modify(~ bind_rows(
      tibble(
        benchmark_source = "separately_collected_market_survey_price",
        unit_value_measure = "raw",
        n_item_observations = nrow(.x),
        weighted_correlation = wtd_cor(.x$log_uv_raw, .x$log_market_price, .x$exp_purch),
        weighted_mae = wtd_mae(.x$log_uv_raw, .x$log_market_price, .x$exp_purch),
        weighted_rmse = wtd_rmse(.x$log_uv_raw, .x$log_market_price, .x$exp_purch)
      ),
      tibble(
        benchmark_source = "separately_collected_market_survey_price",
        unit_value_measure = "standardized",
        n_item_observations = nrow(.x),
        weighted_correlation = wtd_cor(.x$log_uv_std, .x$log_market_price, .x$exp_purch),
        weighted_mae = wtd_mae(.x$log_uv_std, .x$log_market_price, .x$exp_purch),
        weighted_rmse = wtd_rmse(.x$log_uv_std, .x$log_market_price, .x$exp_purch)
      )
    )) %>%
    ungroup()

  paired_value_audit <- tibble(
    statistic = c(
      "n_item_rows",
      "share_item_rows_changed_by_paired_winsorization",
      "expenditure_share_changed_by_paired_winsorization",
      "maximum_paired_identity_error",
      "n_households"
    ),
    value = c(
      nrow(item_data),
      mean(item_data$paired_cleaning_changed, na.rm = TRUE),
      sum(item_data$exp_purch[item_data$paired_cleaning_changed], na.rm = TRUE) /
        sum(item_data$exp_purch, na.rm = TRUE),
      max(abs(item_data$paired_identity_error), na.rm = TRUE),
      nrow(household_data)
    )
  )

  household_summary <- household_data %>%
    summarise(
      n_households = n(),
      mean_C_h = mean(C_h, na.rm = TRUE),
      sd_C_h = sd(C_h, na.rm = TRUE),
      mean_abs_C_h_centered = mean(abs(C_h_centered), na.rm = TRUE),
      survey_weighted_mean_C_h = wtd_mean(C_h, hh_wgt),
      survey_weighted_sd_C_h = {
        weighted_mean <- wtd_mean(C_h, hh_wgt)
        weighted_variance <- wtd_mean((C_h - weighted_mean)^2, hh_wgt)
        if (is.finite(weighted_variance)) sqrt(weighted_variance) else NA_real_
      },
      survey_weighted_mean_abs_C_h_centered = wtd_mean(abs_C_h_centered_weighted, hh_wgt),
      p10_C_h_centered = safe_quantile(C_h_centered, 0.10),
      p50_C_h_centered = safe_quantile(C_h_centered, 0.50),
      p90_C_h_centered = safe_quantile(C_h_centered, 0.90),
      mean_nonkg_share_matched = mean(share_nonkg_exp_matched, na.rm = TRUE),
      median_nonkg_share_matched = median(share_nonkg_exp_matched, na.rm = TRUE),
      mean_conversion_coverage = mean(matched_exp_share_clean, na.rm = TRUE),
      median_conversion_coverage = median(matched_exp_share_clean, na.rm = TRUE)
    )

  item_measurement_construction_manifest <- tribble(
    ~step_order, ~output_variable, ~construction_rule, ~cleaning_or_fallback, ~weighting,
    1L, "log_uv_std_original", "log(exp_purch / (q_purch_raw * factor))", "Positive expenditure, quantity, and conversion factor required", "None",
    2L, "log_uv_std", "Item-level common-unit log value", paste0("Winsorized within item_code at p=", config$item_winsor_p, " in each tail"), "None",
    3L, "log_uv_raw", "log_uv_std + log(factor)", "Reconstructed after common-unit cleaning; not independently winsorized", "None",
    4L, "C_h", "sum(exp_purch * log(factor)) / sum(exp_purch)", "Conversion-matched purchased foods", "Purchase expenditure",
    5L, "share_nonkg_exp_matched", "Expenditure in original unit codes other than kilogram divided by matched expenditure", "Original survey unit_code_original is used", "Purchase expenditure",
    6L, "market benchmark", "Compare item log unit value with separately collected kilogram-based Market Survey log price", "Finite matched pairs only", "Purchase expenditure"
  )

  saveRDS(
    item_data,
    file.path(config$empirical_dir, "analysis_item_level_paired_uv.rds")
  )
  saveRDS(
    household_data,
    file.path(config$empirical_dir, "household_measurement_unit.rds")
  )
  if (isTRUE(config$write_local_microdata_csv)) {
    readr::write_csv(
      household_data,
      file.path(config$empirical_dir, "household_measurement_unit.csv")
    )
  }
  readr::write_csv(
    item_dictionary,
    file.path(config$empirical_dir, "food_item_dictionary.csv")
  )
  readr::write_csv(
    item_winsorization_thresholds,
    file.path(config$empirical_dir, "item_level_winsorization_thresholds.csv")
  )
  readr::write_csv(
    item_measurement_construction_manifest,
    file.path(config$empirical_dir, "item_measurement_construction_manifest.csv")
  )
  readr::write_csv(
    factor_source_summary,
    file.path(config$empirical_dir, "conversion_factor_sources_detailed.csv")
  )
  readr::write_csv(
    conversion_rule_examples,
    file.path(config$empirical_dir, "conversion_rule_examples.csv")
  )
  readr::write_csv(
    market_match_all,
    file.path(config$empirical_dir, "market_price_match_all_conversion_matched.csv")
  )
  readr::write_csv(
    unit_value_market_benchmark,
    file.path(config$empirical_dir, "unit_value_market_price_benchmark.csv")
  )
  readr::write_csv(
    benchmark_by_match_level,
    file.path(config$empirical_dir, "unit_value_market_price_benchmark_by_match_level.csv")
  )
  readr::write_csv(
    household_summary,
    file.path(config$empirical_dir, "household_measurement_unit_summary.csv")
  )
  readr::write_csv(
    paired_value_audit,
    file.path(config$audit_dir, "paired_unit_value_and_C_h_audit.csv")
  )

  list(
    item_data = item_data,
    household_data = household_data,
    item_dictionary = item_dictionary,
    unit_value_market_benchmark = unit_value_market_benchmark,
    item_winsorization_thresholds = item_winsorization_thresholds,
    item_measurement_manifest = item_measurement_construction_manifest
  )
}

# Empirical estimation
estimate_empirical <- function(config, prepared) {
  item_data <- prepared$item_data
  household_data <- prepared$household_data
  item_dictionary <- prepared$item_dictionary

  make_safe_names <- function(x) {
    x %>%
      stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
      stringr::str_replace_all("_+", "_") %>%
      stringr::str_replace_all("^_|_$", "") %>%
      stringr::str_to_lower()
  }

  wtd_sd <- function(x, w) {
    weighted_mean <- wtd_mean(x, w)
    weighted_variance <- wtd_mean((x - weighted_mean)^2, w)
    if (!is.finite(weighted_variance)) NA_real_ else sqrt(weighted_variance)
  }

  summarise_index_difference <- function(
      difference,
      weights,
      all_food_C_centered = NULL,
      system_C_centered = NULL,
      index_type = NA_character_,
      scenario = NULL) {

    difference <- as.numeric(difference)
    weights <- ifelse(is.finite(weights) & weights > 0, weights, 1)
    unweighted_mean <- mean(difference, na.rm = TRUE)
    weighted_mean <- wtd_mean(difference, weights)
    centered_unweighted <- difference - unweighted_mean
    centered_weighted <- difference - weighted_mean

    output <- tibble(
      index_type = index_type,
      n_households = sum(is.finite(difference)),
      mean_difference = unweighted_mean,
      survey_weighted_mean_difference = weighted_mean,
      sd_difference = stats::sd(difference, na.rm = TRUE),
      survey_weighted_sd_difference = wtd_sd(difference, weights),
      mean_absolute_difference_uncentered = mean(abs(difference), na.rm = TRUE),
      survey_weighted_mean_absolute_difference_uncentered =
        wtd_mean(abs(difference), weights),
      mean_absolute_centered_difference = mean(abs(centered_unweighted), na.rm = TRUE),
      survey_weighted_mean_absolute_centered_difference =
        wtd_mean(abs(centered_weighted), weights),
      p10_centered_difference = safe_quantile(centered_unweighted, 0.10),
      p50_centered_difference = safe_quantile(centered_unweighted, 0.50),
      p90_centered_difference = safe_quantile(centered_unweighted, 0.90),
      correlation_with_all_food_C_h_centered = if (is.null(all_food_C_centered)) {
        NA_real_
      } else {
        suppressWarnings(cor(difference, all_food_C_centered, use = "complete.obs"))
      },
      survey_weighted_correlation_with_all_food_C_h_centered =
        if (is.null(all_food_C_centered)) {
          NA_real_
        } else {
          wtd_cor(difference, all_food_C_centered, weights)
        },
      correlation_with_system_C_h_centered = if (is.null(system_C_centered)) {
        NA_real_
      } else {
        suppressWarnings(cor(difference, system_C_centered, use = "complete.obs"))
      },
      survey_weighted_correlation_with_system_C_h_centered =
        if (is.null(system_C_centered)) {
          NA_real_
        } else {
          wtd_cor(difference, system_C_centered, weights)
        }
    )

    if (!is.null(scenario)) output <- mutate(output, scenario = scenario, .before = 1)
    output
  }

  matrix_from_wide <- function(data, prefix, group_names_safe) {
    columns <- paste0(prefix, group_names_safe)
    missing_columns <- setdiff(columns, names(data))
    if (length(missing_columns) > 0) {
      stop(
        "Missing wide columns for prefix ", prefix, ": ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    result <- as.matrix(data[, columns, drop = FALSE])
    storage.mode(result) <- "double"
    colnames(result) <- group_names_safe
    result
  }

  load_group_mapping <- function() {
    mapping <- readr::read_csv(config$group_mapping_path, show_col_types = FALSE) %>%
      mutate(item_code = as.integer(item_code))

    if (!"food_group" %in% names(mapping)) {
      stop("The food-group mapping does not contain a food_group column.", call. = FALSE)
    }

    required <- c("item_code", "item_name", "food_group", "include_main6", "include_main8", "include_all_clean")
    if (!all(required %in% names(mapping))) {
      stop("The food-group mapping is missing required columns.", call. = FALSE)
    }
    mapping <- mapping %>% mutate(across(starts_with("include_"), as.logical))

    duplicated_item_codes <- mapping %>%
      group_by(item_code) %>%
      summarise(n_rows = n(), n_groups = n_distinct(food_group), .groups = "drop") %>%
      filter(n_rows > 1)
    if (nrow(duplicated_item_codes) > 0) {
      stop(
        "The food-group mapping contains duplicated item_code values. Each item must map ",
        "to exactly one row before joining to purchases. See item code(s): ",
        paste(head(duplicated_item_codes$item_code, 20), collapse = ", "),
        call. = FALSE
      )
    }

    mapping
  }

  select_group_mapping <- function(mapping, group_set = "main6", drop_groups = character(0)) {
    include_column <- paste0("include_", group_set)
    if (!include_column %in% names(mapping)) {
      stop("Unknown group set or missing mapping column: ", include_column, call. = FALSE)
    }

    mapping %>%
      mutate(
        include_selected = as.logical(.data[[include_column]]),
        food_group = as.character(food_group)
      ) %>%
      filter(
        include_selected,
        !is.na(food_group), food_group != "",
        !(food_group %in% drop_groups)
      ) %>%
      select(item_code, item_name, food_group)
  }

  read_market_price_files <- function() {
    rm_path <- file.path(config$input_cache_dir, "market_price_item_region_month.csv")
    nm_path <- file.path(config$input_cache_dir, "market_price_item_national_month.csv")
    na_path <- file.path(config$input_cache_dir, "market_price_item_national_allmonths.csv")
    invisible(lapply(c(rm_path, nm_path, na_path), stop_if_missing))

    market_prices <- list(
      region_month = readr::read_csv(rm_path, show_col_types = FALSE) %>%
        transmute(
          item_code = as.integer(item_code),
          region = as.integer(mrk_region),
          food_ym = as.Date(market_ym),
          log_market_rm = as.numeric(log_market_price_rm)
        ),
      national_month = readr::read_csv(nm_path, show_col_types = FALSE) %>%
        transmute(
          item_code = as.integer(item_code),
          food_ym = as.Date(market_ym),
          log_market_nm = as.numeric(log_market_price_nm)
        ),
      national_allmonths = readr::read_csv(na_path, show_col_types = FALSE) %>%
        transmute(
          item_code = as.integer(item_code),
          log_market_na = as.numeric(log_market_price_n)
        )
    )

    assert_unique_key(
      market_prices$region_month,
      c("item_code", "region", "food_ym"),
      "market_price_item_region_month"
    )
    assert_unique_key(
      market_prices$national_month,
      c("item_code", "food_ym"),
      "market_price_item_national_month"
    )
    assert_unique_key(
      market_prices$national_allmonths,
      "item_code",
      "market_price_item_national_allmonths"
    )

    market_prices
  }

  fill_price_hierarchy <- function(data, variable, prefix) {
    variable_name <- rlang::as_name(rlang::ensym(variable))

    med_rm <- data %>%
      group_by(food_group, region, food_ym) %>%
      summarise(!!paste0(prefix, "_rm") := safe_median(.data[[variable_name]]), .groups = "drop")
    med_r <- data %>%
      group_by(food_group, region) %>%
      summarise(!!paste0(prefix, "_r") := safe_median(.data[[variable_name]]), .groups = "drop")
    med_nm <- data %>%
      group_by(food_group, food_ym) %>%
      summarise(!!paste0(prefix, "_nm") := safe_median(.data[[variable_name]]), .groups = "drop")
    med_n <- data %>%
      group_by(food_group) %>%
      summarise(!!paste0(prefix, "_n") := safe_median(.data[[variable_name]]), .groups = "drop")

    generated <- c(
      paste0(prefix, "_rm"), paste0(prefix, "_r"),
      paste0(prefix, "_nm"), paste0(prefix, "_n")
    )
    filled_name <- paste0(variable_name, "_filled")
    fill_level_name <- paste0(variable_name, "_fill_level")

    data %>%
      left_join(med_rm, by = c("food_group", "region", "food_ym")) %>%
      left_join(med_r, by = c("food_group", "region")) %>%
      left_join(med_nm, by = c("food_group", "food_ym")) %>%
      left_join(med_n, by = "food_group") %>%
      mutate(
        !!fill_level_name := case_when(
          is.finite(.data[[variable_name]]) ~ "observed_household_group",
          is.finite(.data[[paste0(prefix, "_rm")]]) ~ "region_month_median",
          is.finite(.data[[paste0(prefix, "_r")]]) ~ "region_all_months_median",
          is.finite(.data[[paste0(prefix, "_nm")]]) ~ "national_month_median",
          is.finite(.data[[paste0(prefix, "_n")]]) ~ "national_all_months_median",
          TRUE ~ "unfilled"
        ),
        !!filled_name := coalesce(
          .data[[variable_name]],
          .data[[paste0(prefix, "_rm")]],
          .data[[paste0(prefix, "_r")]],
          .data[[paste0(prefix, "_nm")]],
          .data[[paste0(prefix, "_n")]]
        )
      ) %>%
      select(-all_of(generated))
  }

  make_price_index <- function(L, W, type = "current", weights = NULL) {
    type <- match.arg(type, c("current", "corrected", "base", "tornqvist"))
    N <- nrow(L)
    G <- ncol(L)
    if (is.null(weights)) weights <- rep(1, N)
    weights <- ifelse(is.finite(weights) & weights > 0, weights, 1)

    base_w <- colSums(sweep(W, 1, weights, `*`), na.rm = TRUE) / sum(weights)
    base_w <- base_w / sum(base_w)
    base_l <- vapply(seq_len(G), function(j) wtd_mean(L[, j], weights), numeric(1))
    base_w_matrix <- matrix(base_w, nrow = N, ncol = G, byrow = TRUE)
    base_l_matrix <- matrix(base_l, nrow = N, ncol = G, byrow = TRUE)

    index <- switch(
      type,
      current = rowSums(W * L),
      corrected = rowSums(W * (L - base_l_matrix)),
      base = as.vector((L - base_l_matrix) %*% base_w),
      tornqvist = rowSums(0.5 * (W + base_w_matrix) * (L - base_l_matrix))
    )

    list(index = index, base_w = base_w, base_l = base_l)
  }

  index_derivative_components <- function(index_type, L_eval, W_eval, base_w, base_l) {
    index_type <- match.arg(index_type, c("current", "corrected", "base", "tornqvist"))
    if (index_type == "current") {
      return(list(direct = W_eval, cvec = W_eval * L_eval))
    }
    if (index_type == "corrected") {
      L_difference <- L_eval - base_l
      return(list(direct = W_eval, cvec = W_eval * L_difference))
    }
    if (index_type == "base") {
      return(list(direct = base_w, cvec = rep(0, length(W_eval))))
    }
    L_difference <- L_eval - base_l
    list(
      direct = 0.5 * (W_eval + base_w),
      cvec = 0.5 * W_eval * L_difference
    )
  }

  build_group_dataset <- function(
      item_data,
      household_data,
      group_mapping,
      scenario_name = "baseline_main6_cov90",
      group_set = "main6",
      min_conversion_coverage = 0.90,
      market_mode = c("fallback", "region_month_only"),
      allowed_factor_sources = NULL,
      excluded_factor_sources = character(0),
      drop_groups = character(0),
      group_winsor_p = config$group_winsor_p) {

    market_mode <- match.arg(market_mode)
    selected_mapping <- select_group_mapping(group_mapping, group_set, drop_groups)
    if (nrow(selected_mapping) == 0) stop("No items selected for scenario ", scenario_name, call. = FALSE)

    require_columns(
      item_data,
      c(
        "case_id", "item_code", "region", "food_ym", "exp_purch",
        "log_uv_std", "log_factor", "factor_source_detailed",
        "is_nonkg_reported_unit"
      ),
      "item_data"
    )
    require_columns(
      household_data,
      c(
        "case_id", "region", "food_ym", "food_month", "hhsize",
        "total_food_purchase_all", "urban",
        "head_female", "head_age", "dependency_ratio"
      ),
      "household_data"
    )
    assert_unique_key(household_data, "case_id", "household_data")
    if (!"hh_wgt" %in% names(household_data)) {
      warning(
        "household_data has no hh_wgt column; equal weights will be used. ",
        "Record this fallback in the empirical construction audit.",
        call. = FALSE
      )
      household_data$hh_wgt <- 1
    }

    item_use <- item_data %>%
      mutate(
        case_id = as.character(case_id),
        item_code = as.integer(item_code),
        region = as.integer(region),
        food_ym = as.Date(food_ym)
      )

    if (!is.null(allowed_factor_sources)) {
      item_use <- item_use %>% filter(factor_source_detailed %in% allowed_factor_sources)
    }
    if (length(excluded_factor_sources) > 0) {
      item_use <- item_use %>% filter(!(factor_source_detailed %in% excluded_factor_sources))
    }

    household_data <- household_data %>%
      mutate(
        case_id = as.character(case_id),
        region = as.integer(region),
        food_ym = as.Date(food_ym),
        food_month = as.integer(food_month),
        urban = as.numeric(urban),
        head_female = as.numeric(head_female),
        head_age = as.numeric(head_age),
        dependency_ratio = as.numeric(dependency_ratio),
        hh_wgt = as.numeric(hh_wgt),
        log_hhsize = log(pmax(as.numeric(hhsize), 1))
      )

    scenario_exposure <- item_use %>%
      group_by(case_id) %>%
      summarise(
        scenario_matched_expenditure = sum(exp_purch, na.rm = TRUE),
        C_h_all = sum(exp_purch * log_factor, na.rm = TRUE) / scenario_matched_expenditure,
        share_nonkg_exp_matched = sum(exp_purch[is_nonkg_reported_unit], na.rm = TRUE) /
          scenario_matched_expenditure,
        .groups = "drop"
      )

    household_coverage_sample <- household_data %>%
      select(-any_of(c("C_h_all", "share_nonkg_exp_matched"))) %>%
      left_join(scenario_exposure, by = "case_id") %>%
      mutate(
        scenario_conversion_coverage = scenario_matched_expenditure / total_food_purchase_all
      ) %>%
      filter(
        is.finite(scenario_matched_expenditure), scenario_matched_expenditure > 0,
        is.finite(scenario_conversion_coverage),
        scenario_conversion_coverage >= min_conversion_coverage
      )

    model_input_variables <- c(
      "region", "food_ym", "food_month", "urban", "head_female",
      "head_age", "dependency_ratio", "log_hhsize", "hh_wgt"
    )
    household_sample <- household_coverage_sample %>%
      filter(
        if_all(all_of(model_input_variables), ~ !is.na(.x)),
        if_all(
          all_of(setdiff(model_input_variables, "food_ym")),
          ~ is.finite(as.numeric(.x))
        ),
        hh_wgt > 0
      )

    if (nrow(household_sample) == 0) {
      stop(
        "No households remain after applying conversion coverage, complete-control, ",
        "date, and positive survey-weight requirements in scenario ", scenario_name, ".",
        call. = FALSE
      )
    }

    household_model_input_sample <- household_sample

    household_sample_audit <- tibble(
      scenario = scenario_name,
      stage = c(
        "coverage_eligible_before_model_inputs",
        "complete_model_inputs_and_positive_weight"
      ),
      n_households = c(nrow(household_coverage_sample), nrow(household_sample))
    )

    item_selected <- item_use %>%
      inner_join(selected_mapping %>% select(item_code, food_group), by = "item_code") %>%
      semi_join(household_sample %>% select(case_id), by = "case_id") %>%
      filter(
        is.finite(exp_purch), exp_purch > 0,
        is.finite(log_uv_std), is.finite(log_factor)
      )

    if (nrow(item_selected) == 0) stop("No usable item observations in scenario ", scenario_name, call. = FALSE)

    group_positive <- item_selected %>%
      group_by(case_id, food_group) %>%
      summarise(
        exp_group = sum(exp_purch, na.rm = TRUE),
        log_uv_std_hh = wtd_mean(log_uv_std, exp_purch),
        log_factor_hh = wtd_mean(log_factor, exp_purch),
        log_uv_raw_hh = log_uv_std_hh + log_factor_hh,
        n_items_group = n_distinct(item_code),
        .groups = "drop"
      ) %>%
      left_join(household_sample %>% select(case_id, region, food_ym), by = "case_id")

    base_item_weights <- item_selected %>%
      group_by(food_group, item_code) %>%
      summarise(exp_item_total = sum(exp_purch, na.rm = TRUE), .groups = "drop") %>%
      group_by(food_group) %>%
      mutate(base_item_weight = exp_item_total / sum(exp_item_total, na.rm = TRUE)) %>%
      ungroup()

    market_prices <- read_market_price_files()
    household_item_grid <- tidyr::expand_grid(
      case_id = household_sample$case_id,
      item_code = unique(base_item_weights$item_code)
    ) %>%
      left_join(household_sample %>% select(case_id, region, food_ym), by = "case_id") %>%
      left_join(base_item_weights, by = "item_code") %>%
      left_join(market_prices$region_month, by = c("item_code", "region", "food_ym"))

    if (market_mode == "fallback") {
      household_item_grid <- household_item_grid %>%
        left_join(market_prices$national_month, by = c("item_code", "food_ym")) %>%
        left_join(market_prices$national_allmonths, by = "item_code") %>%
        mutate(
          log_market_item = coalesce(log_market_rm, log_market_nm, log_market_na),
          market_item_match_level = case_when(
            is.finite(log_market_rm) ~ "region_month",
            !is.finite(log_market_rm) & is.finite(log_market_nm) ~ "national_month",
            !is.finite(log_market_rm) & !is.finite(log_market_nm) & is.finite(log_market_na) ~ "national_allmonths",
            TRUE ~ "unmatched"
          )
        )
    } else {
      household_item_grid <- household_item_grid %>%
        mutate(
          log_market_item = log_market_rm,
          market_item_match_level = if_else(is.finite(log_market_rm), "region_month", "unmatched")
        )
    }

    market_group <- household_item_grid %>%
      group_by(case_id, food_group) %>%
      summarise(
        available_weight = sum(base_item_weight[is.finite(log_market_item)], na.rm = TRUE),
        log_market_base_group = if_else(
          available_weight > 0,
          sum(base_item_weight * log_market_item, na.rm = TRUE) / available_weight,
          NA_real_
        ),
        market_item_weight_coverage = available_weight,
        .groups = "drop"
      )

    groups <- sort(unique(selected_mapping$food_group))
    groups_safe <- make_safe_names(groups)
    if (any(!nzchar(groups_safe)) || anyDuplicated(groups_safe)) {
      stop(
        "Food-group labels do not map one-to-one to nonempty safe column names in scenario ",
        scenario_name, ". Rename the conflicting food_group values.",
        call. = FALSE
      )
    }
    group_name_map <- tibble(food_group = groups, group_safe = groups_safe)

    group_grid <- tidyr::expand_grid(
      case_id = household_sample$case_id,
      food_group = groups
    ) %>%
      left_join(household_sample %>% select(case_id, region, food_ym), by = "case_id") %>%
      left_join(group_positive, by = c("case_id", "food_group", "region", "food_ym")) %>%
      left_join(market_group, by = c("case_id", "food_group")) %>%
      mutate(
        exp_group = replace_na(exp_group, 0),
        n_items_group = replace_na(n_items_group, 0L)
      )

    group_full <- group_grid %>%
      fill_price_hierarchy(log_uv_std_hh, "uvstd") %>%
      fill_price_hierarchy(log_factor_hh, "logfactor")

    if (market_mode == "fallback") {
      group_full <- group_full %>%
        fill_price_hierarchy(log_market_base_group, "market") %>%
        mutate(
          log_market_group_filled = log_market_base_group_filled,
          log_market_group_fill_level = log_market_base_group_fill_level
        )
    } else {
      group_full <- group_full %>%
        mutate(
          log_market_group_filled = log_market_base_group,
          log_market_group_fill_level = if_else(
            is.finite(log_market_base_group),
            "item_aggregation_region_month_only",
            "unfilled"
          )
        )
    }

    price_complete_households <- group_full %>%
      group_by(case_id) %>%
      summarise(
        n_group_rows = n(),
        total_selected_group_expenditure = sum(exp_group, na.rm = TRUE),
        all_group_inputs_finite = all(
          is.finite(log_uv_std_hh_filled) &
            is.finite(log_factor_hh_filled) &
            is.finite(log_market_group_filled)
        ),
        .groups = "drop"
      ) %>%
      filter(
        n_group_rows == length(groups),
        is.finite(total_selected_group_expenditure),
        total_selected_group_expenditure > 0,
        all_group_inputs_finite
      ) %>%
      select(case_id)

    group_full <- group_full %>%
      semi_join(price_complete_households, by = "case_id")
    household_sample <- household_sample %>%
      semi_join(price_complete_households, by = "case_id")

    if (nrow(household_sample) == 0) {
      stop(
        "No household has a complete set of group unit values, conversion factors, and market prices in scenario ",
        scenario_name, ".",
        call. = FALSE
      )
    }

    group_winsorization_thresholds <- bind_rows(
      winsor_thresholds_by(
        group_full,
        log_uv_std_hh_filled,
        by = food_group,
        p = group_winsor_p,
        stage = "group_common_unit_log_value"
      ),
      winsor_thresholds_by(
        group_full,
        log_market_group_filled,
        by = food_group,
        p = group_winsor_p,
        stage = "group_market_survey_log_price"
      )
    ) %>%
      mutate(scenario = scenario_name, .before = 1)

    group_full <- group_full %>%
      winsorize_by(log_uv_std_hh_filled, by = food_group, p = group_winsor_p) %>%
      rename(log_uv_std_group_clean = log_uv_std_hh_filled) %>%
      mutate(
        log_factor_group = log_factor_hh_filled,
        log_uv_raw_group_clean = log_uv_std_group_clean + log_factor_group
      ) %>%
      winsorize_by(log_market_group_filled, by = food_group, p = group_winsor_p) %>%
      group_by(case_id) %>%
      mutate(
        total_group_exp = sum(exp_group, na.rm = TRUE),
        w_group = if_else(total_group_exp > 0, exp_group / total_group_exp, NA_real_),
        n_groups_positive = sum(exp_group > 0, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      left_join(group_name_map, by = "food_group") %>%
      left_join(household_sample %>% select(-region, -food_ym), by = "case_id") %>%
      filter(is.finite(total_group_exp), total_group_exp > 0)

    paired_group_error <- with(
      group_full,
      (log_uv_raw_group_clean - log_uv_std_group_clean) - log_factor_group
    )
    max_group_error <- max(abs(paired_group_error), na.rm = TRUE)
    if (!is.finite(max_group_error) || max_group_error > 1e-10) {
      stop("Group-level paired identity failed in scenario ", scenario_name, call. = FALSE)
    }

    shares_wide <- group_full %>%
      select(case_id, group_safe, w_group) %>%
      pivot_wider(names_from = group_safe, values_from = w_group, values_fill = 0, names_prefix = "w_")
    uv_std_wide <- group_full %>%
      select(case_id, group_safe, log_uv_std_group_clean) %>%
      pivot_wider(names_from = group_safe, values_from = log_uv_std_group_clean, names_prefix = "lp_uv_std_")
    uv_raw_wide <- group_full %>%
      select(case_id, group_safe, log_uv_raw_group_clean) %>%
      pivot_wider(names_from = group_safe, values_from = log_uv_raw_group_clean, names_prefix = "lp_uv_raw_")
    factor_wide <- group_full %>%
      select(case_id, group_safe, log_factor_group) %>%
      pivot_wider(names_from = group_safe, values_from = log_factor_group, names_prefix = "lc_factor_")
    market_wide <- group_full %>%
      select(case_id, group_safe, log_market_group_filled) %>%
      pivot_wider(names_from = group_safe, values_from = log_market_group_filled, names_prefix = "lp_market_")

    wide <- household_sample %>%
      left_join(shares_wide, by = "case_id") %>%
      left_join(uv_std_wide, by = "case_id") %>%
      left_join(uv_raw_wide, by = "case_id") %>%
      left_join(factor_wide, by = "case_id") %>%
      left_join(market_wide, by = "case_id") %>%
      left_join(
        group_full %>% distinct(case_id, total_group_exp, n_groups_positive),
        by = "case_id"
      ) %>%
      mutate(log_total_group_exp = log(total_group_exp))

    share_columns <- paste0("w_", groups_safe)
    price_columns <- c(
      paste0("lp_uv_std_", groups_safe),
      paste0("lp_uv_raw_", groups_safe),
      paste0("lc_factor_", groups_safe),
      paste0("lp_market_", groups_safe)
    )

    wide <- wide %>%
      filter(
        is.finite(log_total_group_exp),
        if_all(all_of(share_columns), ~ is.finite(.x)),
        if_all(all_of(price_columns), ~ is.finite(.x))
      )

    household_sample_audit <- bind_rows(
      household_sample_audit,
      tibble(
        scenario = scenario_name,
        stage = "complete_group_shares_and_all_required_price_series",
        n_households = nrow(wide)
      )
    )

    W <- matrix_from_wide(wide, "w_", groups_safe)
    C_group_matrix <- matrix_from_wide(wide, "lc_factor_", groups_safe)
    C_h_system <- rowSums(W * C_group_matrix)

    wide <- wide %>%
      mutate(
        C_h_system = C_h_system,
        C_h_all_centered = C_h_all - mean(C_h_all, na.rm = TRUE),
        C_h_system_centered = C_h_system - mean(C_h_system, na.rm = TRUE),
        scenario = scenario_name
      )

    group_full_final <- group_full %>%
      semi_join(wide %>% select(case_id), by = "case_id")

    item_main_sample <- item_selected %>%
      semi_join(wide %>% select(case_id), by = "case_id")
    market_match_item <- item_main_sample
    if (market_mode == "region_month_only") {
      if (!"log_market_price_rm" %in% names(market_match_item)) {
        stop("log_market_price_rm is unavailable for the region-month-only scenario.", call. = FALSE)
      }
      market_match_item$market_price_match_level <- ifelse(
        is.finite(market_match_item$log_market_price_rm),
        "region_month", "unmatched"
      )
    } else if ("market_price_match_level" %in% names(market_match_item)) {
      market_match_item$market_price_match_level <- as.character(market_match_item$market_price_match_level)
    } else {
      market_match_item$market_price_match_level <- ifelse(
        is.finite(market_match_item$log_market_price),
        "matched_unspecified", "unmatched"
      )
    }
    market_match_summary <- market_match_item %>%
      group_by(market_price_match_level) %>%
      summarise(
        n_rows = n(),
        expenditure = sum(exp_purch, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        row_share = n_rows / sum(n_rows),
        expenditure_share = expenditure / sum(expenditure),
        scenario = scenario_name,
        .before = 1
      )

    group_summary <- group_full_final %>%
      group_by(food_group) %>%
      summarise(
        n_households = n_distinct(case_id),
        n_positive = sum(exp_group > 0, na.rm = TRUE),
        purchase_rate = mean(exp_group > 0, na.rm = TRUE),
        survey_weighted_purchase_rate = wtd_mean(as.numeric(exp_group > 0), hh_wgt),
        total_expenditure = sum(exp_group, na.rm = TRUE),
        mean_share = mean(w_group, na.rm = TRUE),
        unweighted_mean_share = mean(w_group, na.rm = TRUE),
        survey_weighted_mean_share = wtd_mean(w_group, hh_wgt),
        mean_abs_log_factor_positive = mean(abs(log_factor_hh[exp_group > 0]), na.rm = TRUE),
        mean_market_item_weight_coverage = mean(market_item_weight_coverage, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(expenditure_share = total_expenditure / sum(total_expenditure))

    fill_level_summary <- bind_rows(
      group_full_final %>%
        count(food_group, fill_level = log_uv_std_hh_fill_level, name = "n_households") %>%
        mutate(variable = "group_common_unit_log_value"),
      group_full_final %>%
        count(food_group, fill_level = log_factor_hh_fill_level, name = "n_households") %>%
        mutate(variable = "group_log_conversion_factor"),
      group_full_final %>%
        count(food_group, fill_level = log_market_group_fill_level, name = "n_households") %>%
        mutate(variable = "group_market_survey_log_price")
    ) %>%
      group_by(variable, food_group) %>%
      mutate(household_share = n_households / sum(n_households)) %>%
      ungroup() %>%
      mutate(scenario = scenario_name, .before = 1)

    construction_manifest <- tribble(
      ~step_order, ~output_object, ~construction_rule, ~fallback_or_cleaning, ~weighting,
      1L, "group expenditure", "Sum item purchase expenditure within household x food group", "Zero for non-purchased groups", "None",
      2L, "group standardized log unit value", "sum(item expenditure x item standardized log unit value) / group expenditure", "For non-purchasers: region-month, region-all-months, national-month, national-all-months medians", "Item purchase expenditure",
      3L, "group log conversion factor", "sum(item expenditure x log conversion factor) / group expenditure", "Same hierarchical median fill for non-purchasers", "Item purchase expenditure",
      4L, "group raw log unit value", "cleaned group standardized log unit value + filled group log conversion factor", "Reconstructed after common-unit winsorization", "None",
      5L, "base item weights for market price", "Pooled item expenditure in the coverage-eligible, complete-control-input scenario sample divided by pooled group expenditure", "Fixed before the final price-completeness filter and renormalized over items with an available market price", "Pooled item purchase expenditure",
      6L, "item market price hierarchy", "Region-month, then national-month, then national-all-months", paste0("market_mode=", market_mode), "None",
      7L, "group market log price", "Weighted mean of available item market log prices", "Available item weights renormalized; group-level hierarchical median fill in fallback mode", "Fixed base item weights",
      8L, "group-level winsorization", "Common-unit group log unit value and market-survey group log price", paste0("p=", group_winsor_p, " in each tail within food group; thresholds are estimated on the pre-final household-group grid actually used to construct prices"), "None",
      9L, "budget share", "group purchase expenditure divided by total purchase expenditure in included groups", "Zero shares retained", "None"
    ) %>%
      mutate(scenario = scenario_name, .before = 1)

    list(
      scenario = scenario_name,
      group_set = group_set,
      selected_mapping = selected_mapping,
      groups = groups,
      groups_safe = groups_safe,
      household_coverage_sample = household_coverage_sample,
      household_model_input_sample = household_model_input_sample,
      household_sample = household_sample,
      household_sample_audit = household_sample_audit,
      item_selected = item_main_sample,
      group_long = group_full_final,
      wide = wide,
      group_summary = group_summary,
      market_match_summary = market_match_summary,
      fill_level_summary = fill_level_summary,
      group_winsorization_thresholds = group_winsorization_thresholds,
      base_item_weights = base_item_weights,
      construction_manifest = construction_manifest,
      max_group_identity_error = max_group_error
    )
  }

  price_elasticity_partial <- function(gamma, wbar) {
    groups <- names(wbar)
    output <- expand_grid(i = groups, j = groups) %>%
      mutate(
        elasticity = map2_dbl(i, j, function(ii, jj) {
          -as.numeric(ii == jj) + gamma[ii, jj] / wbar[ii]
        })
      )
    output
  }

  price_elasticity_LA_prime <- function(gamma, beta, wbar, direct) {
    groups <- names(wbar)
    expand_grid(i = groups, j = groups) %>%
      mutate(
        elasticity = map2_dbl(i, j, function(ii, jj) {
          -as.numeric(ii == jj) + gamma[ii, jj] / wbar[ii] -
            (beta[ii] / wbar[ii]) * direct[jj]
        })
      )
  }

  price_elasticity_GA <- function(gamma, beta, wbar, direct, cvec) {
    groups <- names(wbar)
    n <- length(groups)
    bvec <- beta[groups] / wbar[groups]
    A <- diag(n) + as.numeric(bvec) %o% as.numeric(cvec[groups])
    dimnames(A) <- list(groups, groups)

    result <- matrix(NA_real_, n, n, dimnames = list(groups, groups))
    for (j in groups) {
      a <- vapply(groups, function(i) {
        -as.numeric(i == j) + gamma[i, j] / wbar[i] -
          bvec[i] * (direct[j] + cvec[j])
      }, numeric(1))
      solution <- try(solve(A, a), silent = TRUE)
      if (!inherits(solution, "try-error")) result[, j] <- as.numeric(solution)
    }

    grid <- expand.grid(
      i = rownames(result),
      j = colnames(result),
      stringsAsFactors = FALSE
    )
    grid$elasticity <- as.vector(result)
    tibble::as_tibble(grid)
  }

  expenditure_elasticities <- function(beta, wbar, cvec) {
    groups <- names(wbar)
    bvec <- beta[groups] / wbar[groups]
    conventional <- 1 + bvec

    A <- diag(length(groups)) + as.numeric(bvec) %o% as.numeric(cvec[groups])
    adjusted_z <- try(solve(A, as.numeric(bvec)), silent = TRUE)
    adjusted <- if (inherits(adjusted_z, "try-error")) {
      rep(NA_real_, length(groups))
    } else {
      1 + as.numeric(adjusted_z)
    }

    bind_rows(
      tibble(
        group = groups,
        expenditure_elasticity_definition = "conventional_1_plus_beta_over_share",
        expenditure_elasticity = as.numeric(conventional)
      ),
      tibble(
        group = groups,
        expenditure_elasticity_definition = "Green_Alston_current_weight_adjusted",
        expenditure_elasticity = adjusted
      )
    )
  }

  prepare_laaids_sample <- function(
      data,
      groups_safe,
      price_sources,
      use_controls = TRUE,
      require_unique_case_id = FALSE,
      context = "LA/AIDS analysis") {

    data_work <- data
    if (!"hh_wgt" %in% names(data_work)) {
      warning(
        context, " has no hh_wgt column; equal weights will be used.",
        call. = FALSE
      )
      data_work$hh_wgt <- 1
    }

    price_prefix <- c(
      market = "lp_market_",
      uv_std = "lp_uv_std_",
      uv_raw = "lp_uv_raw_"
    )
    unknown_sources <- setdiff(price_sources, names(price_prefix))
    if (length(unknown_sources) > 0) {
      stop(
        "Unknown price source(s) in ", context, ": ",
        paste(unknown_sources, collapse = ", "),
        call. = FALSE
      )
    }

    control_variables <- if (isTRUE(use_controls)) {
      c(
        "urban", "log_hhsize", "head_female", "head_age",
        "dependency_ratio", "region", "food_month"
      )
    } else {
      character(0)
    }
    required_columns <- unique(c(
      "case_id", "log_total_group_exp", "hh_wgt",
      paste0("w_", groups_safe),
      unlist(lapply(price_sources, function(source) {
        paste0(unname(price_prefix[source]), groups_safe)
      }), use.names = FALSE),
      control_variables
    ))
    require_columns(data_work, required_columns, context)

    if (isTRUE(require_unique_case_id)) {
      assert_unique_key(data_work, "case_id", context)
    }

    complete_rows <- stats::complete.cases(data_work[, required_columns, drop = FALSE])
    numeric_columns <- required_columns[
      vapply(data_work[, required_columns, drop = FALSE], is.numeric, logical(1))
    ]
    for (column in numeric_columns) {
      complete_rows <- complete_rows & is.finite(data_work[[column]])
    }
    complete_rows <- complete_rows & data_work$hh_wgt > 0

    analysis_data <- data_work[complete_rows, , drop = FALSE]
    if (nrow(analysis_data) == 0) {
      stop(
        "No complete observations remain for ", context,
        " after applying the common sample and positive-weight requirements.",
        call. = FALSE
      )
    }

    W_check <- matrix_from_wide(analysis_data, "w_", groups_safe)
    maximum_share_sum_error <- max(abs(rowSums(W_check) - 1), na.rm = TRUE)
    if (!is.finite(maximum_share_sum_error) || maximum_share_sum_error > 1e-8) {
      stop(
        "Budget shares do not sum to one in ", context,
        ". Maximum absolute row-sum error = ",
        format(maximum_share_sum_error, scientific = TRUE),
        call. = FALSE
      )
    }

    sample_audit <- tibble(
      context = context,
      n_input_rows = nrow(data_work),
      n_analysis_rows = nrow(analysis_data),
      n_rows_dropped_for_common_sample = nrow(data_work) - nrow(analysis_data),
      common_sample_share = nrow(analysis_data) / nrow(data_work),
      use_controls = isTRUE(use_controls),
      price_sources_required = paste(unique(price_sources), collapse = ";"),
      maximum_share_sum_error = maximum_share_sum_error
    )

    list(
      data = analysis_data,
      audit = sample_audit,
      control_variables = control_variables
    )
  }

  fit_laaids <- function(
      data,
      groups_safe,
      spec_name,
      price_reg = c("market", "uv_std", "uv_raw"),
      price_index = c("market", "uv_std", "uv_raw"),
      index_type = c("current", "corrected", "base", "tornqvist"),
      omitted_group = NULL,
      use_controls = TRUE) {

    price_reg <- match.arg(price_reg)
    price_index <- match.arg(price_index)
    index_type <- match.arg(index_type)

    prefix_reg <- switch(
      price_reg,
      market = "lp_market_",
      uv_std = "lp_uv_std_",
      uv_raw = "lp_uv_raw_"
    )
    prefix_index <- switch(
      price_index,
      market = "lp_market_",
      uv_std = "lp_uv_std_",
      uv_raw = "lp_uv_raw_"
    )

    prepared_sample <- prepare_laaids_sample(
      data = data,
      groups_safe = groups_safe,
      price_sources = unique(c(price_reg, price_index)),
      use_controls = use_controls,
      require_unique_case_id = FALSE,
      context = paste0("specification ", spec_name)
    )
    analysis_data <- prepared_sample$data

    W <- matrix_from_wide(analysis_data, "w_", groups_safe)
    L_reg <- matrix_from_wide(analysis_data, prefix_reg, groups_safe)
    L_index <- matrix_from_wide(analysis_data, prefix_index, groups_safe)
    weights <- as.numeric(analysis_data$hh_wgt)

    index_details <- make_price_index(L_index, W, type = index_type, weights = weights)
    data2 <- analysis_data %>% mutate(
      real_expenditure = log_total_group_exp - index_details$index,
      hh_wgt_safe = weights
    )

    if (is.null(omitted_group)) omitted_group <- groups_safe[length(groups_safe)]
    omitted_group <- make_safe_names(omitted_group)
    if (!omitted_group %in% groups_safe) {
      stop("Omitted group not found: ", omitted_group, call. = FALSE)
    }
    estimated_groups <- setdiff(groups_safe, omitted_group)

    for (group in estimated_groups) {
      data2[[paste0("relative_price_", group)]] <-
        L_reg[, group] - L_reg[, omitted_group]
    }
    relative_price_terms <- paste0("relative_price_", estimated_groups)

    control_variables <- prepared_sample$control_variables
    control_terms <- character(0)
    if (isTRUE(use_controls)) {
      control_terms <- c(
        "urban", "log_hhsize", "head_female", "head_age", "dependency_ratio"
      )
      if (dplyr::n_distinct(data2$region) > 1L) {
        control_terms <- c(control_terms, "factor(region)")
      } else {
        warning(
          "Region fixed effects are constant and therefore omitted in ", spec_name, ".",
          call. = FALSE
        )
      }
      if (dplyr::n_distinct(data2$food_month) > 1L) {
        control_terms <- c(control_terms, "factor(food_month)")
      } else {
        warning(
          "Interview-month fixed effects are constant and therefore omitted in ",
          spec_name, ".",
          call. = FALSE
        )
      }
    }

    coefficient_rows <- list()
    model_stat_rows <- list()
    beta <- setNames(rep(NA_real_, length(groups_safe)), groups_safe)
    gamma <- matrix(
      0,
      nrow = length(groups_safe), ncol = length(groups_safe),
      dimnames = list(groups_safe, groups_safe)
    )

    for (group in estimated_groups) {
      dependent <- paste0("w_", group)
      rhs <- c(relative_price_terms, "real_expenditure", control_terms)
      formula <- stats::as.formula(paste(dependent, "~", paste(rhs, collapse = " + ")))

      selected_columns <- unique(c(
        dependent, relative_price_terms, "real_expenditure", "hh_wgt_safe",
        control_variables
      ))
      require_columns(data2, selected_columns, paste0(spec_name, " equation ", group))
      model_data <- data2 %>% select(all_of(selected_columns))
      if (any(!stats::complete.cases(model_data))) {
        stop(
          "The common-sample preparation left missing values in ", spec_name,
          ", equation ", group, ".",
          call. = FALSE
        )
      }
      numeric_model_columns <- names(model_data)[
        vapply(model_data, is.numeric, logical(1))
      ]
      if (any(vapply(
        model_data[numeric_model_columns],
        function(x) any(!is.finite(x)),
        logical(1)
      ))) {
        stop(
          "The common-sample preparation left non-finite numeric values in ",
          spec_name, ", equation ", group, ".",
          call. = FALSE
        )
      }

      fit <- stats::lm(formula, data = model_data, weights = hh_wgt_safe)
      if (stats::nobs(fit) != nrow(data2)) {
        stop(
          "Equation-specific row deletion occurred in ", spec_name,
          ", equation ", group, ". All equations must use the common sample.",
          call. = FALSE
        )
      }
      coefficient_rows[[group]] <- broom::tidy(fit) %>%
        mutate(spec = spec_name, equation = group, .before = 1)
      model_stat_rows[[group]] <- broom::glance(fit) %>%
        mutate(spec = spec_name, equation = group, n_obs_equation = stats::nobs(fit), .before = 1)

      coefficients <- stats::coef(fit)
      beta_value <- unname(coefficients["real_expenditure"])
      if (!is.finite(beta_value)) {
        stop("Non-finite real-expenditure coefficient in ", spec_name, ", equation ", group, call. = FALSE)
      }
      beta[group] <- beta_value
      for (j in estimated_groups) {
        term <- paste0("relative_price_", j)
        coefficient_value <- if (term %in% names(coefficients)) {
          unname(coefficients[term])
        } else {
          NA_real_
        }
        if (!is.finite(coefficient_value)) {
          stop(
            "Non-finite relative-price coefficient ", term, " in ",
            spec_name, ", equation ", group, call. = FALSE
          )
        }
        gamma[group, j] <- coefficient_value
      }
      gamma[group, omitted_group] <- -sum(gamma[group, estimated_groups])
    }

    beta[omitted_group] <- -sum(beta[estimated_groups], na.rm = TRUE)
    for (j in groups_safe) {
      gamma[omitted_group, j] <- -sum(gamma[estimated_groups, j], na.rm = TRUE)
    }

    W_weighted <- sweep(W, 1, weights, `*`)
    wbar <- colSums(W_weighted, na.rm = TRUE) / sum(weights)
    wbar <- wbar / sum(wbar)
    names(wbar) <- groups_safe

    L_index_eval <- vapply(
      seq_along(groups_safe),
      function(j) wtd_mean(L_index[, j], weights),
      numeric(1)
    )
    names(L_index_eval) <- groups_safe
    names(index_details$base_w) <- groups_safe
    names(index_details$base_l) <- groups_safe

    derivative <- index_derivative_components(
      index_type = index_type,
      L_eval = L_index_eval,
      W_eval = wbar,
      base_w = index_details$base_w,
      base_l = index_details$base_l
    )
    names(derivative$direct) <- groups_safe
    names(derivative$cvec) <- groups_safe

    expenditure_elasticities <- expenditure_elasticities(
      beta = beta,
      wbar = wbar,
      cvec = derivative$cvec
    ) %>%
      mutate(
        spec = spec_name,
        beta = beta[group],
        mean_share = wbar[group],
        price_regressor = price_reg,
        index_price = price_index,
        index_type = index_type,
        .before = 1
      )

    partial <- price_elasticity_partial(gamma, wbar) %>%
      mutate(elasticity_definition = "partial_regressor_index_fixed")
    LA_prime <- price_elasticity_LA_prime(
      gamma, beta, wbar, derivative$direct
    ) %>%
      mutate(elasticity_definition = "joint_one_for_one_LA_prime")
    Green_Alston <- price_elasticity_GA(
      gamma, beta, wbar, derivative$direct, derivative$cvec
    ) %>%
      mutate(elasticity_definition = "joint_one_for_one_Green_Alston")

    price_elasticities <- bind_rows(partial, LA_prime, Green_Alston) %>%
      mutate(
        spec = spec_name,
        elasticity_type = if_else(i == j, "own", "cross"),
        gamma_ij = map2_dbl(i, j, ~ gamma[.x, .y]),
        beta_i = beta[i],
        w_i = wbar[i],
        w_j = wbar[j],
        price_regressor = price_reg,
        index_price = price_index,
        index_type = index_type,
        .before = 1
      )

    symmetry_audit <- tibble(
      spec = spec_name,
      n_input_rows = prepared_sample$audit$n_input_rows,
      n_analysis_rows = prepared_sample$audit$n_analysis_rows,
      n_rows_dropped_for_common_sample =
        prepared_sample$audit$n_rows_dropped_for_common_sample,
      maximum_share_sum_error = prepared_sample$audit$maximum_share_sum_error,
      maximum_absolute_gamma_asymmetry = max(abs(gamma - t(gamma)), na.rm = TRUE),
      mean_absolute_gamma_asymmetry = mean(abs(gamma - t(gamma)), na.rm = TRUE),
      adding_up_beta_error = abs(sum(beta)),
      maximum_gamma_row_sum_error = max(abs(rowSums(gamma))),
      maximum_gamma_column_sum_error = max(abs(colSums(gamma)))
    )

    list(
      spec = spec_name,
      coefficients = bind_rows(coefficient_rows),
      model_stats = bind_rows(model_stat_rows),
      beta = beta,
      gamma = gamma,
      wbar = wbar,
      expenditure_elasticities = expenditure_elasticities,
      price_elasticities = price_elasticities,
      symmetry_audit = symmetry_audit,
      sample_audit = prepared_sample$audit,
      omitted_group = omitted_group,
      data_n = nrow(data2),
      index_details = index_details,
      derivative = derivative
    )
  }

  compare_expenditure_specs <- function(expenditure_data, raw_spec, std_spec, comparison) {
    expenditure_data %>%
      filter(spec %in% c(raw_spec, std_spec)) %>%
      select(
        spec, group, expenditure_elasticity_definition,
        beta, mean_share, expenditure_elasticity
      ) %>%
      pivot_wider(
        names_from = spec,
        values_from = c(beta, expenditure_elasticity)
      ) %>%
      mutate(
        comparison = comparison,
        beta_raw_minus_std = .data[[paste0("beta_", raw_spec)]] -
          .data[[paste0("beta_", std_spec)]],
        expenditure_elasticity_raw_minus_std =
          .data[[paste0("expenditure_elasticity_", raw_spec)]] -
          .data[[paste0("expenditure_elasticity_", std_spec)]]
      )
  }

  compare_price_specs <- function(price_data, raw_spec, std_spec, comparison) {
    price_data %>%
      filter(spec %in% c(raw_spec, std_spec)) %>%
      select(spec, i, j, elasticity_type, elasticity_definition, elasticity) %>%
      pivot_wider(names_from = spec, values_from = elasticity) %>%
      mutate(
        comparison = comparison,
        raw_minus_std = .data[[raw_spec]] - .data[[std_spec]],
        absolute_difference = abs(raw_minus_std),
        sign_change = sign(.data[[raw_spec]]) != sign(.data[[std_spec]])
      )
  }

  detect_bootstrap_variable <- function(
      data, requested, candidates, role, allow_none = FALSE) {
    requested <- if (is.null(requested)) "" else trimws(as.character(requested)[1])
    if (nzchar(requested)) {
      if (!requested %in% names(data)) {
        stop(
          role, " variable requested through the environment was not found: ",
          requested,
          call. = FALSE
        )
      }
      return(requested)
    }

    available <- candidates[candidates %in% names(data)]
    if (length(available) > 0) return(available[1])
    if (allow_none) return(NA_character_)
    if ("case_id" %in% names(data)) return("case_id")
    stop("No usable ", role, " variable was found.", call. = FALSE)
  }

  prepare_empirical_bootstrap_design <- function(data, cfg) {
    require_columns(data, "case_id", "bootstrap analysis data")
    assert_unique_key(data, "case_id", "bootstrap analysis data")

    cluster_variable <- detect_bootstrap_variable(
      data = data,
      requested = cfg$empirical_bootstrap_cluster_var,
      candidates = c(
        "ea_id", "ea", "cluster_id", "cluster", "psu", "psu_id",
        "enumeration_area", "enumeration_area_id"
      ),
      role = "bootstrap cluster",
      allow_none = FALSE
    )
    strata_variable <- detect_bootstrap_variable(
      data = data,
      requested = cfg$empirical_bootstrap_strata_var,
      candidates = c(
        "strata", "stratum", "strata_id", "stratum_id",

        "district"
      ),
      role = "bootstrap stratum",
      allow_none = TRUE
    )

    if (is.na(strata_variable) &&
        isTRUE(cfg$empirical_bootstrap_require_strata)) {
      stop(
        "No survey-stratum variable was found for the empirical bootstrap. ",
        "The Malawi IHS5 used 32 strata. Ensure that `district` is retained in ",
        "hh_base_controls.csv and check the district field. ",
        "The paper requires a stratified ",
        "diagnostic run, not for the final manuscript.",
        call. = FALSE
      )
    }

    bootstrap_data <- data
    cluster_values <- as.character(bootstrap_data[[cluster_variable]])
    missing_cluster <- is.na(cluster_values) | trimws(cluster_values) == ""
    cluster_values[missing_cluster] <- paste0(
      "__singleton_missing_cluster__", bootstrap_data$case_id[missing_cluster]
    )
    cluster_is_household_unique <- n_distinct(cluster_values) == nrow(bootstrap_data)
    effective_household_level_fallback <- identical(cluster_variable, "case_id") ||
      all(missing_cluster) || cluster_is_household_unique

    if (is.na(strata_variable)) {
      strata_values <- rep("__all_households__", nrow(bootstrap_data))
    } else {
      strata_values <- as.character(bootstrap_data[[strata_variable]])
      strata_values[is.na(strata_values) | trimws(strata_values) == ""] <-
        "__missing_stratum__"
    }

    bootstrap_data$.boot_stratum <- strata_values
    bootstrap_data$.boot_cluster_original <- cluster_values
    bootstrap_data$.boot_cluster_key <- paste(strata_values, cluster_values, sep = "::")

    cluster_sizes <- bootstrap_data %>%
      count(.boot_stratum, .boot_cluster_key, name = "n_households_cluster")
    clusters_per_stratum <- cluster_sizes %>%
      count(.boot_stratum, name = "n_clusters_stratum")
    minimum_clusters_per_stratum <- min(clusters_per_stratum$n_clusters_stratum)
    observed_strata <- n_distinct(bootstrap_data$.boot_stratum)
    expected_strata <- cfg$empirical_bootstrap_expected_strata
    strata_count_matches_expected <- is.finite(expected_strata) &&
      observed_strata == expected_strata

    if (isTRUE(cfg$empirical_bootstrap_require_strata) && observed_strata < 2L) {
      stop(
        "The final empirical bootstrap resolved to only ", observed_strata,
        " stratum. The Malawi IHS5 final analysis should be resampled within ",
        "its district sampling strata. Inspect district/stratum propagation ",
        "before using bootstrap intervals.",
        call. = FALSE
      )
    }

    if (is.finite(expected_strata) && observed_strata != expected_strata) {
      warning(
        "The empirical bootstrap contains ", observed_strata,
        " observed strata, whereas expected strata = ", expected_strata,
        ". This may be legitimate if an analysis restriction removes an entire ",
        "stratum, but it must be checked before reporting the intervals.",
        call. = FALSE
      )
    }

    if (effective_household_level_fallback) {
      warning(
        "No usable non-missing PSU/EA variable was detected for the empirical bootstrap. ",
        "Household-level resampling will be used. Check the survey ",
        "PSU variable when it is available in hh_base_controls.csv.",
        call. = FALSE
      )
    } else if (any(missing_cluster)) {
      warning(
        sum(missing_cluster), " household(s) have a missing bootstrap cluster value and ",
        "will be treated as singleton clusters. Inspect empirical_bootstrap_design.csv.",
        call. = FALSE
      )
    }
    if (!is.na(strata_variable) && minimum_clusters_per_stratum < 2L) {
      warning(
        "At least one bootstrap stratum contains fewer than two clusters. Such strata ",
        "cannot contribute resampling variation; inspect empirical_bootstrap_design.csv.",
        call. = FALSE
      )
    }

    design <- tibble(
      bootstrap_method = "paired nonparametric cluster bootstrap within strata",
      uncertainty_scope = paste(
        "Conditional on the constructed household-group analysis file;",
        "the same resampled households are used for raw and standardized specifications;",
        "price-index reference means/share weights (where applicable) and every model are re-estimated;",
        "item/group winsorization thresholds, conversion factors, hierarchical fills, and constructed household variables are held fixed."
      ),
      cluster_variable = cluster_variable,
      strata_variable = if (is.na(strata_variable)) "none" else strata_variable,
      strata_required = isTRUE(cfg$empirical_bootstrap_require_strata),
      expected_number_of_strata = expected_strata,
      observed_number_of_strata = observed_strata,
      observed_strata_match_expected = strata_count_matches_expected,
      household_level_fallback = effective_household_level_fallback,
      selected_cluster_variable_is_household_unique = cluster_is_household_unique,
      n_missing_original_cluster_values = sum(missing_cluster),
      share_missing_original_cluster_values = mean(missing_cluster),
      n_households = nrow(bootstrap_data),
      n_clusters = nrow(cluster_sizes),
      n_strata = n_distinct(bootstrap_data$.boot_stratum),
      minimum_clusters_per_stratum = minimum_clusters_per_stratum,
      median_clusters_per_stratum = stats::median(clusters_per_stratum$n_clusters_stratum),
      maximum_clusters_per_stratum = max(clusters_per_stratum$n_clusters_stratum),
      minimum_cluster_size = min(cluster_sizes$n_households_cluster),
      median_cluster_size = stats::median(cluster_sizes$n_households_cluster),
      maximum_cluster_size = max(cluster_sizes$n_households_cluster),
      requested_replications = cfg$empirical_bootstrap_reps,
      seed = cfg$empirical_bootstrap_seed,
      percentile_interval = "2.5th and 97.5th percentiles",
      minimum_required_success_rate = cfg$empirical_bootstrap_min_success_rate
    )

    list(data = bootstrap_data, design = design)
  }

  resample_empirical_clusters <- function(bootstrap_data) {
    cluster_table <- bootstrap_data %>%
      distinct(.boot_stratum, .boot_cluster_key)
    clusters_by_stratum <- split(cluster_table, cluster_table$.boot_stratum, drop = TRUE)

    draw_counter <- 0L
    sampled_rows <- lapply(clusters_by_stratum, function(stratum_clusters) {
      sampled_index <- sample.int(
        n = nrow(stratum_clusters),
        size = nrow(stratum_clusters),
        replace = TRUE
      )
      bind_rows(lapply(sampled_index, function(index) {
        draw_counter <<- draw_counter + 1L
        key <- stratum_clusters$.boot_cluster_key[index]
        rows <- bootstrap_data[
          bootstrap_data$.boot_cluster_key == key,
          ,
          drop = FALSE
        ]
        rows$.boot_cluster_draw <- draw_counter
        rows
      }))
    })

    bind_rows(sampled_rows)
  }

  run_one_empirical_bootstrap <- function(
      bootstrap_data, groups_safe, omitted_group, replication, replication_seed) {
    set.seed(replication_seed)
    sampled <- resample_empirical_clusters(bootstrap_data)

    W_sampled <- matrix_from_wide(sampled, "w_", groups_safe)
    L_uv_raw_sampled <- matrix_from_wide(sampled, "lp_uv_raw_", groups_safe)
    L_uv_std_sampled <- matrix_from_wide(sampled, "lp_uv_std_", groups_safe)
    bootstrap_weights <- ifelse(
      is.finite(sampled$hh_wgt) & sampled$hh_wgt > 0,
      sampled$hh_wgt,
      1
    )
    index_sensitivity <- map_dfr(
      c("current", "corrected", "tornqvist", "base"),
      function(index_type) {
        raw_details <- make_price_index(
          L_uv_raw_sampled, W_sampled,
          type = index_type, weights = bootstrap_weights
        )
        standardized_details <- make_price_index(
          L_uv_std_sampled, W_sampled,
          type = index_type, weights = bootstrap_weights
        )
        difference <- raw_details$index - standardized_details$index
        summarise_index_difference(
          difference = difference,
          weights = bootstrap_weights,
          all_food_C_centered = sampled$C_h_all,
          system_C_centered = sampled$C_h_system,
          index_type = index_type
        ) %>%
          mutate(replication = replication, .before = 1)
      }
    )

    central_specs <- tribble(
      ~spec,                     ~price_reg, ~price_index, ~index_type,
      "LA_PM_UVstd_current",     "market",   "uv_std",     "current",
      "LA_PM_UVraw_current",     "market",   "uv_raw",     "current"
    )
    fits <- pmap(
      central_specs,
      function(spec, price_reg, price_index, index_type) {
        fit_laaids(
          data = sampled,
          groups_safe = groups_safe,
          spec_name = spec,
          price_reg = price_reg,
          price_index = price_index,
          index_type = index_type,
          omitted_group = omitted_group,
          use_controls = TRUE
        )
      }
    )

    expenditure_levels <- map_dfr(fits, "expenditure_elasticities") %>%
      select(
        spec, group, expenditure_elasticity_definition,
        beta, mean_share, expenditure_elasticity
      ) %>%
      mutate(replication = replication, .before = 1)

    expenditure_comparison <- compare_expenditure_specs(
      map_dfr(fits, "expenditure_elasticities"),
      raw_spec = "LA_PM_UVraw_current",
      std_spec = "LA_PM_UVstd_current",
      comparison = "current_stone_only_raw_vs_standardized"
    ) %>%
      transmute(
        replication = replication,
        group,
        expenditure_elasticity_definition,
        beta_raw = .data[["beta_LA_PM_UVraw_current"]],
        beta_standardized = .data[["beta_LA_PM_UVstd_current"]],
        beta_raw_minus_standardized = beta_raw_minus_std,
        expenditure_elasticity_raw =
          .data[["expenditure_elasticity_LA_PM_UVraw_current"]],
        expenditure_elasticity_standardized =
          .data[["expenditure_elasticity_LA_PM_UVstd_current"]],
        expenditure_elasticity_raw_minus_standardized =
          expenditure_elasticity_raw_minus_std,
        opposite_sides_of_one =
          (expenditure_elasticity_raw - 1) *
          (expenditure_elasticity_standardized - 1) < 0
      )

    price_levels <- map_dfr(fits, "price_elasticities") %>%
      select(
        spec, i, j, elasticity_type, elasticity_definition,
        elasticity, gamma_ij, beta_i, w_i, w_j
      ) %>%
      mutate(replication = replication, .before = 1)

    price_comparison <- compare_price_specs(
      map_dfr(fits, "price_elasticities"),
      raw_spec = "LA_PM_UVraw_current",
      std_spec = "LA_PM_UVstd_current",
      comparison = "current_stone_only_raw_vs_standardized"
    ) %>%
      transmute(
        replication = replication,
        i, j, elasticity_type, elasticity_definition,
        elasticity_raw = .data[["LA_PM_UVraw_current"]],
        elasticity_standardized = .data[["LA_PM_UVstd_current"]],
        elasticity_raw_minus_standardized = raw_minus_std,
        sign_change = sign_change
      )

    finite_checks <- c(
      index_sensitivity$mean_difference,
      index_sensitivity$sd_difference,
      index_sensitivity$mean_absolute_centered_difference,
      index_sensitivity$survey_weighted_sd_difference,
      index_sensitivity$survey_weighted_mean_absolute_centered_difference,
      expenditure_levels$beta,
      expenditure_levels$mean_share,
      expenditure_levels$expenditure_elasticity,
      expenditure_comparison$beta_raw,
      expenditure_comparison$beta_standardized,
      expenditure_comparison$beta_raw_minus_standardized,
      expenditure_comparison$expenditure_elasticity_raw,
      expenditure_comparison$expenditure_elasticity_standardized,
      expenditure_comparison$expenditure_elasticity_raw_minus_standardized,
      price_levels$elasticity,
      price_comparison$elasticity_raw,
      price_comparison$elasticity_standardized,
      price_comparison$elasticity_raw_minus_standardized
    )
    if (any(!is.finite(finite_checks))) {
      stop(
        "A central bootstrap fit returned a non-finite coefficient or elasticity; ",
        "the replication is recorded as failed rather than silently entering an interval.",
        call. = FALSE
      )
    }

    list(
      index_sensitivity = index_sensitivity,
      expenditure_levels = expenditure_levels,
      expenditure_comparison = expenditure_comparison,
      price_levels = price_levels,
      price_comparison = price_comparison,
      n_resampled_household_rows = nrow(sampled),
      n_cluster_draws = n_distinct(sampled$.boot_cluster_draw),
      n_unique_original_clusters_represented = n_distinct(sampled$.boot_cluster_key)
    )
  }

  run_empirical_bootstrap <- function(data, groups_safe, omitted_group, cfg) {
    prepared <- prepare_empirical_bootstrap_design(data, cfg)
    bootstrap_data <- prepared$data
    n_replications <- cfg$empirical_bootstrap_reps

    index_sensitivity_rows <- vector("list", n_replications)
    expenditure_level_rows <- vector("list", n_replications)
    expenditure_comparison_rows <- vector("list", n_replications)
    price_level_rows <- vector("list", n_replications)
    price_comparison_rows <- vector("list", n_replications)
    diagnostic_rows <- vector("list", n_replications)

    for (replication in seq_len(n_replications)) {
      if (replication == 1L ||
          replication %% cfg$empirical_bootstrap_progress_every == 0L ||
          replication == n_replications) {
        message(
          "Empirical paired bootstrap replication ", replication,
          " of ", n_replications
        )
      }

      replication_seed <- as.integer(cfg$empirical_bootstrap_seed + replication - 1L)
      result <- tryCatch(
        run_one_empirical_bootstrap(
          bootstrap_data = bootstrap_data,
          groups_safe = groups_safe,
          omitted_group = omitted_group,
          replication = replication,
          replication_seed = replication_seed
        ),
        error = function(e) e
      )

      if (inherits(result, "error")) {
        diagnostic_rows[[replication]] <- tibble(
          replication = replication,
          replication_seed = replication_seed,
          status = "failed",
          n_resampled_household_rows = NA_integer_,
          n_cluster_draws = NA_integer_,
          n_unique_original_clusters_represented = NA_integer_,
          error_message = conditionMessage(result)
        )
        next
      }

      index_sensitivity_rows[[replication]] <- result$index_sensitivity
      expenditure_level_rows[[replication]] <- result$expenditure_levels
      expenditure_comparison_rows[[replication]] <- result$expenditure_comparison
      price_level_rows[[replication]] <- result$price_levels
      price_comparison_rows[[replication]] <- result$price_comparison
      diagnostic_rows[[replication]] <- tibble(
        replication = replication,
        replication_seed = replication_seed,
        status = "completed",
        n_resampled_household_rows = result$n_resampled_household_rows,
        n_cluster_draws = result$n_cluster_draws,
        n_unique_original_clusters_represented =
          result$n_unique_original_clusters_represented,
        error_message = NA_character_
      )
    }

    diagnostics <- bind_rows(diagnostic_rows)
    successful_replications <- sum(diagnostics$status == "completed")
    success_rate <- successful_replications / n_replications
    design <- prepared$design %>%
      mutate(
        successful_replications = successful_replications,
        failed_replications = n_replications - successful_replications,
        success_rate = success_rate
      )

    list(
      design = design,
      diagnostics = diagnostics,
      index_sensitivity = bind_rows(index_sensitivity_rows),
      expenditure_levels = bind_rows(expenditure_level_rows),
      expenditure_comparison = bind_rows(expenditure_comparison_rows),
      price_levels = bind_rows(price_level_rows),
      price_comparison = bind_rows(price_comparison_rows)
    )
  }

  summarise_empirical_bootstrap <- function(
      bootstrap, point_index_sensitivity,
      point_expenditure_levels, point_expenditure_comparison,
      point_price_levels, point_price_comparison) {

    index_metrics <- c(
      "mean_difference",
      "survey_weighted_mean_difference",
      "sd_difference",
      "survey_weighted_sd_difference",
      "mean_absolute_difference_uncentered",
      "survey_weighted_mean_absolute_difference_uncentered",
      "mean_absolute_centered_difference",
      "survey_weighted_mean_absolute_centered_difference",
      "correlation_with_all_food_C_h_centered",
      "survey_weighted_correlation_with_all_food_C_h_centered",
      "correlation_with_system_C_h_centered",
      "survey_weighted_correlation_with_system_C_h_centered"
    )
    index_point_long <- point_index_sensitivity %>%
      select(index_type, all_of(index_metrics)) %>%
      pivot_longer(
        cols = all_of(index_metrics),
        names_to = "metric",
        values_to = "point_estimate"
      )
    index_sensitivity_summary <- bootstrap$index_sensitivity %>%
      select(replication, index_type, all_of(index_metrics)) %>%
      pivot_longer(
        cols = all_of(index_metrics),
        names_to = "metric",
        values_to = "replication_statistic"
      ) %>%
      group_by(index_type, metric) %>%
      summarise(
        n_successful_bootstrap = n_distinct(replication),
        bootstrap_mean = mean(replication_statistic, na.rm = TRUE),
        bootstrap_se = sd(replication_statistic, na.rm = TRUE),
        ci_lower = safe_quantile(replication_statistic, 0.025),
        ci_upper = safe_quantile(replication_statistic, 0.975),
        .groups = "drop"
      ) %>%
      left_join(index_point_long, by = c("index_type", "metric")) %>%
      mutate(
        bootstrap_method = "paired cluster percentile bootstrap",
        interval = "2.5th and 97.5th percentiles",
        location_or_dispersion = case_when(
          metric %in% c("mean_difference", "survey_weighted_mean_difference") ~
            "location_shift_absorbed_by_intercept",
          stringr::str_detect(metric, "correlation") ~ "mechanism_correlation",
          stringr::str_detect(metric, "uncentered") ~ "uncentered_diagnostic",
          TRUE ~ "household_dispersion_primary"
        )
      )

    expenditure_level_point <- point_expenditure_levels %>%
      rename(
        point_beta = beta,
        point_mean_share = mean_share,
        point_expenditure_elasticity = expenditure_elasticity
      )
    expenditure_level_summary <- bootstrap$expenditure_levels %>%
      group_by(spec, group, expenditure_elasticity_definition) %>%
      summarise(
        n_successful_bootstrap = n_distinct(replication),
        bootstrap_mean_beta = mean(beta, na.rm = TRUE),
        bootstrap_se_beta = sd(beta, na.rm = TRUE),
        beta_ci_lower = safe_quantile(beta, 0.025),
        beta_ci_upper = safe_quantile(beta, 0.975),
        bootstrap_mean_share = mean(mean_share, na.rm = TRUE),
        bootstrap_mean_expenditure_elasticity = mean(expenditure_elasticity, na.rm = TRUE),
        bootstrap_se_expenditure_elasticity = sd(expenditure_elasticity, na.rm = TRUE),
        expenditure_elasticity_ci_lower = safe_quantile(expenditure_elasticity, 0.025),
        expenditure_elasticity_ci_upper = safe_quantile(expenditure_elasticity, 0.975),
        probability_expenditure_elasticity_above_one =
          mean(expenditure_elasticity > 1, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(
        expenditure_level_point,
        by = c("spec", "group", "expenditure_elasticity_definition")
      )

    expenditure_difference_point <- point_expenditure_comparison %>%
      transmute(
        group,
        expenditure_elasticity_definition,
        point_beta_raw = .data[["beta_LA_PM_UVraw_current"]],
        point_beta_standardized = .data[["beta_LA_PM_UVstd_current"]],
        point_beta_raw_minus_standardized = beta_raw_minus_std,
        point_expenditure_elasticity_raw =
          .data[["expenditure_elasticity_LA_PM_UVraw_current"]],
        point_expenditure_elasticity_standardized =
          .data[["expenditure_elasticity_LA_PM_UVstd_current"]],
        point_expenditure_elasticity_raw_minus_standardized =
          expenditure_elasticity_raw_minus_std,
        point_estimates_opposite_sides_of_one =
          (point_expenditure_elasticity_raw - 1) *
          (point_expenditure_elasticity_standardized - 1) < 0
      )
    expenditure_difference_summary <- bootstrap$expenditure_comparison %>%
      group_by(group, expenditure_elasticity_definition) %>%
      summarise(
        n_successful_bootstrap = n_distinct(replication),
        bootstrap_se_beta_difference = sd(beta_raw_minus_standardized, na.rm = TRUE),
        beta_difference_ci_lower = safe_quantile(beta_raw_minus_standardized, 0.025),
        beta_difference_ci_upper = safe_quantile(beta_raw_minus_standardized, 0.975),
        probability_beta_difference_positive =
          mean(beta_raw_minus_standardized > 0, na.rm = TRUE),
        bootstrap_se_expenditure_elasticity_difference =
          sd(expenditure_elasticity_raw_minus_standardized, na.rm = TRUE),
        expenditure_elasticity_difference_ci_lower =
          safe_quantile(expenditure_elasticity_raw_minus_standardized, 0.025),
        expenditure_elasticity_difference_ci_upper =
          safe_quantile(expenditure_elasticity_raw_minus_standardized, 0.975),
        probability_expenditure_elasticity_difference_positive =
          mean(expenditure_elasticity_raw_minus_standardized > 0, na.rm = TRUE),
        probability_opposite_sides_of_one = mean(opposite_sides_of_one, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(
        expenditure_difference_point,
        by = c("group", "expenditure_elasticity_definition")
      )

    price_level_point <- point_price_levels %>%
      rename(point_elasticity = elasticity)
    price_level_summary <- bootstrap$price_levels %>%
      group_by(spec, i, j, elasticity_type, elasticity_definition) %>%
      summarise(
        n_successful_bootstrap = n_distinct(replication),
        bootstrap_mean_elasticity = mean(elasticity, na.rm = TRUE),
        bootstrap_se_elasticity = sd(elasticity, na.rm = TRUE),
        elasticity_ci_lower = safe_quantile(elasticity, 0.025),
        elasticity_ci_upper = safe_quantile(elasticity, 0.975),
        probability_elasticity_positive = mean(elasticity > 0, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(
        price_level_point,
        by = c("spec", "i", "j", "elasticity_type", "elasticity_definition")
      )

    price_difference_point <- point_price_comparison %>%
      transmute(
        i, j, elasticity_type, elasticity_definition,
        point_elasticity_raw = .data[["LA_PM_UVraw_current"]],
        point_elasticity_standardized = .data[["LA_PM_UVstd_current"]],
        point_elasticity_raw_minus_standardized = raw_minus_std,
        point_sign_change = sign_change
      )
    price_difference_summary <- bootstrap$price_comparison %>%
      group_by(i, j, elasticity_type, elasticity_definition) %>%
      summarise(
        n_successful_bootstrap = n_distinct(replication),
        bootstrap_se_difference = sd(elasticity_raw_minus_standardized, na.rm = TRUE),
        difference_ci_lower = safe_quantile(elasticity_raw_minus_standardized, 0.025),
        difference_ci_upper = safe_quantile(elasticity_raw_minus_standardized, 0.975),
        probability_difference_positive =
          mean(elasticity_raw_minus_standardized > 0, na.rm = TRUE),
        probability_sign_change = mean(sign_change, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(
        price_difference_point,
        by = c("i", "j", "elasticity_type", "elasticity_definition")
      )

    point_expenditure_aggregate <- point_expenditure_comparison %>%
      group_by(expenditure_elasticity_definition) %>%
      summarise(
        mean_absolute_beta_difference = mean(abs(beta_raw_minus_std), na.rm = TRUE),
        maximum_absolute_beta_difference = max(abs(beta_raw_minus_std), na.rm = TRUE),
        mean_absolute_expenditure_elasticity_difference =
          mean(abs(expenditure_elasticity_raw_minus_std), na.rm = TRUE),
        maximum_absolute_expenditure_elasticity_difference =
          max(abs(expenditure_elasticity_raw_minus_std), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_longer(
        cols = -expenditure_elasticity_definition,
        names_to = "metric",
        values_to = "point_estimate"
      )

    expenditure_aggregate_summary <- bootstrap$expenditure_comparison %>%
      group_by(replication, expenditure_elasticity_definition) %>%
      summarise(
        mean_absolute_beta_difference = mean(abs(beta_raw_minus_standardized), na.rm = TRUE),
        maximum_absolute_beta_difference = max(abs(beta_raw_minus_standardized), na.rm = TRUE),
        mean_absolute_expenditure_elasticity_difference =
          mean(abs(expenditure_elasticity_raw_minus_standardized), na.rm = TRUE),
        maximum_absolute_expenditure_elasticity_difference =
          max(abs(expenditure_elasticity_raw_minus_standardized), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_longer(
        cols = -c(replication, expenditure_elasticity_definition),
        names_to = "metric",
        values_to = "replication_statistic"
      ) %>%
      group_by(expenditure_elasticity_definition, metric) %>%
      summarise(
        n_successful_bootstrap = n_distinct(replication),
        bootstrap_mean = mean(replication_statistic, na.rm = TRUE),
        bootstrap_se = sd(replication_statistic, na.rm = TRUE),
        ci_lower = safe_quantile(replication_statistic, 0.025),
        ci_upper = safe_quantile(replication_statistic, 0.975),
        .groups = "drop"
      ) %>%
      left_join(
        point_expenditure_aggregate,
        by = c("expenditure_elasticity_definition", "metric")
      ) %>%
      mutate(
        bootstrap_method = "paired cluster percentile bootstrap",
        interval = "2.5th and 97.5th percentiles"
      )

    point_price_aggregate <- point_price_comparison %>%
      group_by(elasticity_definition, elasticity_type) %>%
      summarise(
        mean_absolute_difference = mean(absolute_difference, na.rm = TRUE),
        median_absolute_difference = median(absolute_difference, na.rm = TRUE),
        maximum_absolute_difference = max(absolute_difference, na.rm = TRUE),
        sign_change_rate = mean(sign_change, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_longer(
        cols = -c(elasticity_definition, elasticity_type),
        names_to = "metric",
        values_to = "point_estimate"
      )

    price_aggregate_summary <- bootstrap$price_comparison %>%
      group_by(replication, elasticity_definition, elasticity_type) %>%
      summarise(
        mean_absolute_difference =
          mean(abs(elasticity_raw_minus_standardized), na.rm = TRUE),
        median_absolute_difference =
          median(abs(elasticity_raw_minus_standardized), na.rm = TRUE),
        maximum_absolute_difference =
          max(abs(elasticity_raw_minus_standardized), na.rm = TRUE),
        sign_change_rate = mean(sign_change, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_longer(
        cols = -c(replication, elasticity_definition, elasticity_type),
        names_to = "metric",
        values_to = "replication_statistic"
      ) %>%
      group_by(elasticity_definition, elasticity_type, metric) %>%
      summarise(
        n_successful_bootstrap = n_distinct(replication),
        bootstrap_mean = mean(replication_statistic, na.rm = TRUE),
        bootstrap_se = sd(replication_statistic, na.rm = TRUE),
        ci_lower = safe_quantile(replication_statistic, 0.025),
        ci_upper = safe_quantile(replication_statistic, 0.975),
        .groups = "drop"
      ) %>%
      left_join(
        point_price_aggregate,
        by = c("elasticity_definition", "elasticity_type", "metric")
      ) %>%
      mutate(
        bootstrap_method = "paired cluster percentile bootstrap",
        interval = "2.5th and 97.5th percentiles"
      )

    list(
      index_sensitivity = index_sensitivity_summary,
      expenditure_levels = expenditure_level_summary,
      expenditure_differences = expenditure_difference_summary,
      price_levels = price_level_summary,
      price_differences = price_difference_summary,
      expenditure_aggregate_sensitivity = expenditure_aggregate_summary,
      price_aggregate_sensitivity = price_aggregate_summary
    )
  }

  group_mapping <- load_group_mapping()
  readr::write_csv(
    group_mapping,
    file.path(config$empirical_dir, "food_group_mapping_used.csv")
  )

  baseline_group_data <- build_group_dataset(
    item_data = item_data,
    household_data = household_data,
    group_mapping = group_mapping,
    scenario_name = "baseline_main6_cov90",
    group_set = "main6",
    min_conversion_coverage = config$baseline_min_conversion_coverage,
    market_mode = "fallback",
    excluded_factor_sources = c("rough_liquid_fallback"),
    group_winsor_p = config$group_winsor_p
  )

  saveRDS(
    baseline_group_data$group_long,
    file.path(config$empirical_dir, "group_long_main6.rds")
  )
  saveRDS(
    baseline_group_data$wide,
    file.path(config$empirical_dir, "group_wide_main6.rds")
  )
  readr::write_csv(
    baseline_group_data$group_summary,
    file.path(config$empirical_dir, "table_group_coverage_main6.csv")
  )
  readr::write_csv(
    baseline_group_data$household_sample_audit,
    file.path(config$empirical_dir, "household_model_input_sample_audit_main6.csv")
  )
  readr::write_csv(
    baseline_group_data$fill_level_summary,
    file.path(config$empirical_dir, "group_price_fill_hierarchy_main6.csv")
  )
  readr::write_csv(
    baseline_group_data$group_winsorization_thresholds,
    file.path(config$empirical_dir, "group_winsorization_thresholds_main6.csv")
  )
  readr::write_csv(
    baseline_group_data$base_item_weights,
    file.path(config$empirical_dir, "market_group_base_item_weights_main6.csv")
  )
  readr::write_csv(
    baseline_group_data$construction_manifest,
    file.path(config$empirical_dir, "group_data_construction_manifest_main6.csv")
  )

  group_selection <- group_mapping %>%
    mutate(
      food_group = as.character(food_group),
      include_main6 = as.logical(include_main6),
      include_main8 = as.logical(include_main8)
    ) %>%
    group_by(food_group) %>%
    summarise(
      n_item_codes = n_distinct(item_code),
      included_main6 = any(include_main6, na.rm = TRUE),
      included_main8 = any(include_main8, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      selection_rationale = case_when(
        included_main6 ~ "Included to balance substantive interpretability, purchase prevalence, conversion coverage, and market-price coverage.",
        food_group == "fruits" ~ "Excluded from the conservative six-group baseline; included in the eight-group robustness system.",
        food_group == "dairy" ~ "Excluded from the conservative six-group baseline; included in the eight-group robustness system.",
        food_group %in% c("prepared_vendor_foods", "prepared_snacks") ~ "Excluded because prepared/vendor foods are not closely comparable with the raw-item Market Survey benchmark.",
        food_group %in% c("sweets_snacks") ~ "Excluded from the conservative core-food system and retained only as an optional category.",
        food_group == "alcohol" ~ "Excluded from the food-demand baseline.",
        food_group %in% c("other_unspecified", "review_needed") ~ "Excluded because item content cannot be mapped consistently.",
        TRUE ~ "Excluded from the conservative six-group baseline; documented for mapping review."
      )
    ) %>%
    left_join(
      baseline_group_data$group_summary %>%
        select(
          food_group, purchase_rate, expenditure_share,
          mean_market_item_weight_coverage
        ),
      by = "food_group"
    ) %>%
    arrange(desc(included_main6), desc(included_main8), food_group)
  readr::write_csv(
    group_selection,
    file.path(config$empirical_dir, "table_food_group_selection.csv")
  )

  readr::write_csv(
    baseline_group_data$market_match_summary,
    file.path(config$empirical_dir, "table_market_price_matching_main6_cov90.csv")
  )
  readr::write_csv(
    baseline_group_data$selected_mapping,
    file.path(config$empirical_dir, "food_group_mapping_selected_main6.csv")
  )

  hh_base_count_path <- file.path(config$input_cache_dir, "hh_base_controls.csv")
  all_household_count <- if (file.exists(hh_base_count_path)) {
    n_distinct(readr::read_csv(hh_base_count_path, show_col_types = FALSE)$case_id)
  } else {
    NA_integer_
  }
  sample_flow <- tibble(
    stage = c(
      "IHS5 households in household control file",
      "Households with at least one conversion-matched purchased food",
      paste0("Households meeting conversion coverage >= ", config$baseline_min_conversion_coverage),
      "Coverage-eligible households with complete model inputs and positive survey weights",
      "Households in final complete six-group diagnostic estimation sample"
    ),
    n_households = c(
      all_household_count,
      n_distinct(household_data$case_id),
      n_distinct(baseline_group_data$household_coverage_sample$case_id),
      n_distinct(baseline_group_data$household_model_input_sample$case_id),
      nrow(baseline_group_data$wide)
    ),
    n_selected_item_rows = c(
      NA_integer_,
      nrow(item_data),
      nrow(
        item_data %>%
          filter(factor_source_detailed != "rough_liquid_fallback") %>%
          semi_join(
            baseline_group_data$household_coverage_sample %>% select(case_id),
            by = "case_id"
          )
      ),
      nrow(
        item_data %>%
          filter(factor_source_detailed != "rough_liquid_fallback") %>%
          semi_join(
            baseline_group_data$household_model_input_sample %>% select(case_id),
            by = "case_id"
          )
      ),
      nrow(baseline_group_data$item_selected)
    )
  )
  readr::write_csv(
    sample_flow,
    file.path(config$empirical_dir, "table_sample_flow.csv")
  )

  writeLines(
    c(
      "Empirical analysis scope for the purchased-food demand analysis",
      "",
      "1. The analytical item file contains purchased foods with positive reported purchased quantity and expenditure.",
      "2. Own-produced and gifted quantities are not assigned constructed market prices and do not enter the unit-value, budget-share, or six-group demand-system calculations.",
      "3. The baseline excludes the universal rough litre-to-kilogram fallback; item-specific liquid-density rules are retained, and inclusion of the rough fallback is checked separately.",
      "4. Market-price regressors use separately collected item-region-month Market Survey prices, followed by national-month and national-all-months fallbacks in the baseline.",
      "5. Item log values are aggregated to household-food-group log values using item purchase-expenditure weights; item market log prices use fixed pooled item-expenditure weights renormalized over available prices.",
      "6. Non-purchaser group unit values are filled by region-month, region-all-months, national-month, and national-all-months medians, in that order.",
      "7. The reported expenditure elasticities are conditional on total purchased expenditure in the included food groups, not unconditional income elasticities for all food consumption.",
      "8. The share equations are a diagnostic homogeneity-and-adding-up implementation with zero shares retained; they are not presented as a preferred censored structural demand system for Malawi.",
      "9. Paired percentile intervals condition on the constructed analysis file: the same cluster-resampled households enter raw and standardized specifications, while base shares, reference log prices, indices, models, and elasticities are re-estimated in every successful replication."
    ),
    con = file.path(config$empirical_dir, "analysis_scope.txt")
  )

  wide_baseline <- baseline_group_data$wide
  groups_safe <- baseline_group_data$groups_safe
  W_baseline <- matrix_from_wide(wide_baseline, "w_", groups_safe)
  L_uv_raw_baseline <- matrix_from_wide(wide_baseline, "lp_uv_raw_", groups_safe)
  L_uv_std_baseline <- matrix_from_wide(wide_baseline, "lp_uv_std_", groups_safe)
  weights_baseline <- ifelse(
    is.finite(wide_baseline$hh_wgt) & wide_baseline$hh_wgt > 0,
    wide_baseline$hh_wgt,
    1
  )

  index_sensitivity <- map_dfr(c("current", "corrected", "tornqvist", "base"), function(index_type) {
    raw_details <- make_price_index(
      L_uv_raw_baseline, W_baseline, type = index_type, weights = weights_baseline
    )
    std_details <- make_price_index(
      L_uv_std_baseline, W_baseline, type = index_type, weights = weights_baseline
    )
    difference <- raw_details$index - std_details$index
    summarise_index_difference(
      difference = difference,
      weights = weights_baseline,
      all_food_C_centered = wide_baseline$C_h_all_centered,
      system_C_centered = wide_baseline$C_h_system_centered,
      index_type = index_type
    ) %>%
      mutate(
        primary_dispersion_metric = "mean_absolute_centered_difference",
        maximum_current_identity_error = if (index_type == "current") {
          max(abs(difference - wide_baseline$C_h_system), na.rm = TRUE)
        } else {
          NA_real_
        }
      )
  })

  readr::write_csv(
    index_sensitivity,
    file.path(config$empirical_dir, "table_index_sensitivity.csv")
  )

  specifications <- tribble(
    ~spec,                          ~price_reg, ~price_index, ~index_type,
    "LA_PM_PM_current",             "market",   "market",     "current",
    "LA_PM_UVstd_current",          "market",   "uv_std",     "current",
    "LA_PM_UVraw_current",          "market",   "uv_raw",     "current",
    "LA_UVstd_UVstd_current",       "uv_std",   "uv_std",     "current",
    "LA_UVraw_UVraw_current",       "uv_raw",   "uv_raw",     "current",
    "LA_PM_UVstd_corrected",        "market",   "uv_std",     "corrected",
    "LA_PM_UVraw_corrected",        "market",   "uv_raw",     "corrected",
    "LA_PM_UVstd_tornqvist",        "market",   "uv_std",     "tornqvist",
    "LA_PM_UVraw_tornqvist",        "market",   "uv_raw",     "tornqvist",
    "LA_PM_UVstd_base",             "market",   "uv_std",     "base",
    "LA_PM_UVraw_base",             "market",   "uv_raw",     "base"
  )

  omitted_group <- groups_safe[length(groups_safe)]

  estimation_specification_manifest <- specifications %>%
    mutate(
      dependent_variables = "Six purchased-food-group budget shares; one equation omitted and recovered by adding-up",
      real_expenditure_regressor = "log(total purchased expenditure in included groups) minus stated price index",
      price_parameterization = paste0("Relative ", price_reg, " log prices; omitted group = ", omitted_group),
      scalar_equation = "w_ih = alpha_i + sum_j gamma_ij log(p_jh/p_Gh) + beta_i[log(x_h)-index_h] + controls + error_ih",
      controls = "Urban, log household size, head sex, head age, dependency ratio, region fixed effects, interview-month fixed effects",
      equation_weighting = "Weighted least squares using positive finite hh_wgt; equal weights are used only when hh_wgt is absent from the source household file",
      maintained_restrictions = "Homogeneity and adding-up; symmetry and curvature not imposed; zero shares retained",
      uncertainty = if_else(
        spec %in% c("LA_PM_UVstd_current", "LA_PM_UVraw_current") &
          isTRUE(config$run_empirical_bootstrap),
        "Paired cluster bootstrap percentile intervals; same resampled households for both specifications",
        if_else(
          spec %in% c("LA_PM_UVstd_current", "LA_PM_UVraw_current"),
          "Bootstrap disabled; point estimates only",
          "Point estimates and sensitivity summaries only"
        )
      ),
      interpretation = "Diagnostic LA/AIDS sensitivity exercise, not a preferred censored structural demand estimator"
    )
  readr::write_csv(
    estimation_specification_manifest,
    file.path(config$empirical_dir, "estimation_specification_manifest.csv")
  )

  model_fits <- pmap(
    specifications,
    function(spec, price_reg, price_index, index_type) {
      fit_laaids(
        data = wide_baseline,
        groups_safe = groups_safe,
        spec_name = spec,
        price_reg = price_reg,
        price_index = price_index,
        index_type = index_type,
        omitted_group = omitted_group,
        use_controls = TRUE
      )
    }
  )

  coefficients_all <- map_dfr(model_fits, "coefficients")
  model_stats_all <- map_dfr(model_fits, "model_stats")
  expenditure_all <- map_dfr(model_fits, "expenditure_elasticities")
  price_all <- map_dfr(model_fits, "price_elasticities")
  system_audit_all <- map_dfr(model_fits, "symmetry_audit")

  readr::write_csv(coefficients_all, file.path(config$empirical_dir, "laaids_coefficients.csv"))
  readr::write_csv(model_stats_all, file.path(config$empirical_dir, "laaids_model_statistics.csv"))
  readr::write_csv(expenditure_all, file.path(config$empirical_dir, "laaids_expenditure_elasticities.csv"))
  readr::write_csv(price_all, file.path(config$empirical_dir, "laaids_price_elasticities_all_definitions.csv"))
  readr::write_csv(system_audit_all, file.path(config$audit_dir, "laaids_system_restriction_audit.csv"))

  central_expenditure_levels <- expenditure_all %>%
    filter(spec %in% c("LA_PM_UVstd_current", "LA_PM_UVraw_current")) %>%
    select(
      spec, group, expenditure_elasticity_definition,
      beta, mean_share, expenditure_elasticity
    )

  central_expenditure_comparison <- compare_expenditure_specs(
    expenditure_all,
    raw_spec = "LA_PM_UVraw_current",
    std_spec = "LA_PM_UVstd_current",
    comparison = "current_stone_only_raw_vs_standardized"
  )

  central_price_levels <- price_all %>%
    filter(spec %in% c("LA_PM_UVstd_current", "LA_PM_UVraw_current")) %>%
    select(
      spec, i, j, elasticity_type, elasticity_definition,
      elasticity, gamma_ij, beta_i, w_i, w_j
    )

  central_price_comparison_cells <- compare_price_specs(
    price_all,
    raw_spec = "LA_PM_UVraw_current",
    std_spec = "LA_PM_UVstd_current",
    comparison = "current_stone_only_raw_vs_standardized"
  )

  central_price_comparison_summary <- central_price_comparison_cells %>%
    group_by(comparison, elasticity_definition, elasticity_type) %>%
    summarise(
      mean_absolute_difference = mean(absolute_difference, na.rm = TRUE),
      median_absolute_difference = median(absolute_difference, na.rm = TRUE),
      maximum_absolute_difference = max(absolute_difference, na.rm = TRUE),
      sign_change_rate = mean(sign_change, na.rm = TRUE),
      .groups = "drop"
    )

  comparison_pairs <- tribble(
    ~comparison,                     ~raw_spec,                    ~std_spec,
    "current_stone_only",            "LA_PM_UVraw_current",        "LA_PM_UVstd_current",
    "corrected_stone",               "LA_PM_UVraw_corrected",      "LA_PM_UVstd_corrected",
    "tornqvist",                     "LA_PM_UVraw_tornqvist",      "LA_PM_UVstd_tornqvist",
    "base_share",                    "LA_PM_UVraw_base",           "LA_PM_UVstd_base",
    "full_raw_vs_standardized",      "LA_UVraw_UVraw_current",     "LA_UVstd_UVstd_current"
  )

  all_expenditure_comparisons <- pmap_dfr(
    comparison_pairs,
    function(comparison, raw_spec, std_spec) {
      compare_expenditure_specs(expenditure_all, raw_spec, std_spec, comparison)
    }
  )
  all_expenditure_comparison_summary <- all_expenditure_comparisons %>%
    group_by(comparison, expenditure_elasticity_definition) %>%
    summarise(
      mean_absolute_beta_difference = mean(abs(beta_raw_minus_std), na.rm = TRUE),
      maximum_absolute_beta_difference = max(abs(beta_raw_minus_std), na.rm = TRUE),
      mean_absolute_expenditure_elasticity_difference =
        mean(abs(expenditure_elasticity_raw_minus_std), na.rm = TRUE),
      maximum_absolute_expenditure_elasticity_difference =
        max(abs(expenditure_elasticity_raw_minus_std), na.rm = TRUE),
      .groups = "drop"
    )

  all_price_comparison_cells <- pmap_dfr(
    comparison_pairs,
    function(comparison, raw_spec, std_spec) {
      compare_price_specs(price_all, raw_spec, std_spec, comparison)
    }
  )
  all_price_comparison_summary <- all_price_comparison_cells %>%
    group_by(comparison, elasticity_definition, elasticity_type) %>%
    summarise(
      mean_absolute_difference = mean(absolute_difference, na.rm = TRUE),
      median_absolute_difference = median(absolute_difference, na.rm = TRUE),
      maximum_absolute_difference = max(absolute_difference, na.rm = TRUE),
      sign_change_rate = mean(sign_change, na.rm = TRUE),
      .groups = "drop"
    )

  readr::write_csv(
    central_expenditure_levels,
    file.path(config$empirical_dir, "table_expenditure_response_levels_current_stone.csv")
  )
  readr::write_csv(
    central_expenditure_comparison,
    file.path(config$empirical_dir, "table_expenditure_response_cells_current_stone.csv")
  )
  readr::write_csv(
    all_expenditure_comparison_summary,
    file.path(config$empirical_dir, "table_expenditure_response_sensitivity.csv")
  )
  readr::write_csv(
    central_price_levels,
    file.path(config$empirical_dir, "table_price_elasticity_levels_current_stone.csv")
  )
  readr::write_csv(
    central_price_comparison_cells,
    file.path(config$empirical_dir, "table_price_elasticity_cells_current_stone.csv")
  )
  readr::write_csv(
    central_price_comparison_summary,
    file.path(config$empirical_dir, "table_price_elasticity_sensitivity_current_stone.csv")
  )
  readr::write_csv(
    all_price_comparison_summary,
    file.path(config$empirical_dir, "table_price_elasticity_sensitivity_all_indices.csv")
  )
  readr::write_csv(
    all_expenditure_comparisons,
    file.path(config$empirical_dir, "table_expenditure_response_cells_all_indices.csv")
  )
  readr::write_csv(
    all_price_comparison_cells,
    file.path(config$empirical_dir, "table_price_elasticity_cells_all_indices.csv")
  )

  empirical_bootstrap <- NULL
  empirical_bootstrap_summary <- NULL
  if (isTRUE(config$run_empirical_bootstrap)) {
    empirical_bootstrap <- run_empirical_bootstrap(
      data = wide_baseline,
      groups_safe = groups_safe,
      omitted_group = omitted_group,
      cfg = config
    )

    readr::write_csv(
      empirical_bootstrap$design,
      file.path(config$empirical_dir, "empirical_bootstrap_design.csv")
    )
    readr::write_csv(
      empirical_bootstrap$diagnostics,
      file.path(config$audit_dir, "empirical_bootstrap_diagnostics.csv")
    )
    saveRDS(
      empirical_bootstrap,
      file.path(config$empirical_dir, "empirical_bootstrap_replicates.rds"),
      compress = TRUE
    )

    bootstrap_success_rate <- empirical_bootstrap$design$success_rate[1]
    if (!is.finite(bootstrap_success_rate) ||
        bootstrap_success_rate < config$empirical_bootstrap_min_success_rate) {
      stop(
        "The empirical paired bootstrap success rate was ",
        format(bootstrap_success_rate, digits = 4),
        ", below required success rate = ",
        config$empirical_bootstrap_min_success_rate,
        ". Inspect empirical_bootstrap_diagnostics.csv.",
        call. = FALSE
      )
    }

    empirical_bootstrap_summary <- summarise_empirical_bootstrap(
      bootstrap = empirical_bootstrap,
      point_index_sensitivity = index_sensitivity,
      point_expenditure_levels = central_expenditure_levels,
      point_expenditure_comparison = central_expenditure_comparison,
      point_price_levels = central_price_levels,
      point_price_comparison = central_price_comparison_cells
    )

    readr::write_csv(
      empirical_bootstrap_summary$index_sensitivity,
      file.path(
        config$empirical_dir,
        "empirical_bootstrap_index_sensitivity_with_ci.csv"
      )
    )
    readr::write_csv(
      empirical_bootstrap_summary$expenditure_levels,
      file.path(
        config$empirical_dir,
        "empirical_bootstrap_expenditure_levels_with_ci.csv"
      )
    )
    readr::write_csv(
      empirical_bootstrap_summary$expenditure_differences,
      file.path(
        config$empirical_dir,
        "empirical_bootstrap_expenditure_differences_with_ci.csv"
      )
    )
    readr::write_csv(
      empirical_bootstrap_summary$price_levels,
      file.path(
        config$empirical_dir,
        "empirical_bootstrap_price_levels_with_ci.csv"
      )
    )
    readr::write_csv(
      empirical_bootstrap_summary$price_differences,
      file.path(
        config$empirical_dir,
        "empirical_bootstrap_price_differences_with_ci.csv"
      )
    )
    readr::write_csv(
      empirical_bootstrap_summary$expenditure_aggregate_sensitivity,
      file.path(
        config$empirical_dir,
        "empirical_bootstrap_expenditure_aggregate_sensitivity_with_ci.csv"
      )
    )
    readr::write_csv(
      empirical_bootstrap_summary$price_aggregate_sensitivity,
      file.path(
        config$empirical_dir,
        "empirical_bootstrap_price_aggregate_sensitivity_with_ci.csv"
      )
    )
  } else {
    readr::write_csv(
      tibble(
        bootstrap_method = "paired nonparametric cluster bootstrap within strata",
        status = "disabled",
        requested_replications = config$empirical_bootstrap_reps,
        seed = config$empirical_bootstrap_seed
      ),
      file.path(config$empirical_dir, "empirical_bootstrap_design.csv")
    )
  }

  benchmark_main_sample <- baseline_group_data$item_selected %>%
    filter(
      is.finite(log_market_price),
      is.finite(log_uv_raw),
      is.finite(log_uv_std),
      is.finite(exp_purch), exp_purch > 0
    )
  benchmark_main_summary <- map_dfr(
    c(raw = "log_uv_raw", standardized = "log_uv_std"),
    function(variable) {
      x <- benchmark_main_sample[[variable]]
      y <- benchmark_main_sample$log_market_price
      w <- benchmark_main_sample$exp_purch
      tibble(
        benchmark_source = "separately_collected_market_survey_price",
        unit_value_measure = if (identical(variable, "log_uv_raw")) "raw" else "standardized",
        n_item_observations = length(x),
        expenditure_weighted_correlation = wtd_cor(x, y, w),
        expenditure_weighted_mae = wtd_mae(x, y, w),
        expenditure_weighted_rmse = wtd_rmse(x, y, w)
      )
    }
  )
  readr::write_csv(
    benchmark_main_summary,
    file.path(config$empirical_dir, "table_unit_value_market_benchmark_main6_cov90.csv")
  )

  baseline_identity_audit <- tibble(
    audit = c(
      "maximum_item_level_paired_identity_error",
      "maximum_group_level_paired_identity_error",
      "maximum_current_index_minus_system_C_h_error",
      "n_households_baseline",
      "n_groups_baseline"
    ),
    value = c(
      max(abs(item_data$paired_identity_error), na.rm = TRUE),
      baseline_group_data$max_group_identity_error,
      index_sensitivity$maximum_current_identity_error[index_sensitivity$index_type == "current"],
      nrow(wide_baseline),
      length(groups_safe)
    )
  )
  readr::write_csv(
    baseline_identity_audit,
    file.path(config$audit_dir, "baseline_empirical_identity_audit.csv")
  )

  empirical_results <- list(
    group_mapping = group_mapping,
    baseline_group_data = baseline_group_data,
    fits = model_fits,
    coefficients = coefficients_all,
    model_stats = model_stats_all,
    expenditure = expenditure_all,
    price = price_all,
    index_sensitivity = index_sensitivity,
    central_expenditure_levels = central_expenditure_levels,
    central_expenditure_comparison = central_expenditure_comparison,
    central_price_levels = central_price_levels,
    central_price_comparison_cells = central_price_comparison_cells,
    central_price_comparison_summary = central_price_comparison_summary,
    all_expenditure_comparisons = all_expenditure_comparisons,
    all_price_comparison_cells = all_price_comparison_cells,
    empirical_bootstrap = empirical_bootstrap,
    empirical_bootstrap_summary = empirical_bootstrap_summary,
    benchmark_main_summary = benchmark_main_summary
  )

  run_robustness_and_figures <- function() {
    robustness_scenarios <- list(
      list(
        scenario_name = "baseline_main6_cov90",
        group_set = "main6",
        min_conversion_coverage = config$baseline_min_conversion_coverage,
        market_mode = "fallback",
        allowed_factor_sources = NULL,
        excluded_factor_sources = c("rough_liquid_fallback"),
        drop_groups = character(0),
        group_winsor_p = config$group_winsor_p
      ),
      list(
        scenario_name = "include_rough_liquid_fallback",
        group_set = "main6",
        min_conversion_coverage = 0.90,
        market_mode = "fallback",
        allowed_factor_sources = NULL,
        excluded_factor_sources = character(0),
        drop_groups = character(0),
        group_winsor_p = config$group_winsor_p
      ),
      list(
        scenario_name = "high_conversion_cov95",
        group_set = "main6",
        min_conversion_coverage = 0.95,
        market_mode = "fallback",
        allowed_factor_sources = NULL,
        excluded_factor_sources = c("rough_liquid_fallback"),
        drop_groups = character(0),
        group_winsor_p = config$group_winsor_p
      ),
      list(
        scenario_name = "broader_conversion_cov80",
        group_set = "main6",
        min_conversion_coverage = 0.80,
        market_mode = "fallback",
        allowed_factor_sources = NULL,
        excluded_factor_sources = c("rough_liquid_fallback"),
        drop_groups = character(0),
        group_winsor_p = config$group_winsor_p
      ),
      list(
        scenario_name = "exclude_manual_liquid_units",
        group_set = "main6",
        min_conversion_coverage = 0.90,
        market_mode = "fallback",
        allowed_factor_sources = NULL,
        excluded_factor_sources = c("manual_item_specific_liquid", "rough_liquid_fallback"),
        drop_groups = character(0),
        group_winsor_p = config$group_winsor_p
      ),
      list(
        scenario_name = "official_or_exact_units_only",
        group_set = "main6",
        min_conversion_coverage = 0.90,
        market_mode = "fallback",
        allowed_factor_sources = c(
          "official_direct", "official_subunit_code", "official_other_text",
          "exact_kilogram", "exact_gram"
        ),
        excluded_factor_sources = character(0),
        drop_groups = character(0),
        group_winsor_p = config$group_winsor_p
      ),
      list(
        scenario_name = "region_month_market_only",
        group_set = "main6",
        min_conversion_coverage = 0.90,
        market_mode = "region_month_only",
        allowed_factor_sources = NULL,
        excluded_factor_sources = c("rough_liquid_fallback"),
        drop_groups = character(0),
        group_winsor_p = config$group_winsor_p
      ),
      list(
        scenario_name = "alternative_group_main8",
        group_set = "main8",
        min_conversion_coverage = 0.90,
        market_mode = "fallback",
        allowed_factor_sources = NULL,
        excluded_factor_sources = c("rough_liquid_fallback"),
        drop_groups = character(0),
        group_winsor_p = config$group_winsor_p
      ),
      list(
        scenario_name = "drop_animal_protein",
        group_set = "main6",
        min_conversion_coverage = 0.90,
        market_mode = "fallback",
        allowed_factor_sources = NULL,
        excluded_factor_sources = c("rough_liquid_fallback"),
        drop_groups = c("animal_protein"),
        group_winsor_p = config$group_winsor_p
      ),
      list(
        scenario_name = "drop_oils_sugar_condiments",
        group_set = "main6",
        min_conversion_coverage = 0.90,
        market_mode = "fallback",
        allowed_factor_sources = NULL,
        excluded_factor_sources = c("rough_liquid_fallback"),
        drop_groups = c("oils_sugar_condiments"),
        group_winsor_p = config$group_winsor_p
      ),
      list(
        scenario_name = "no_group_winsorization",
        group_set = "main6",
        min_conversion_coverage = 0.90,
        market_mode = "fallback",
        allowed_factor_sources = NULL,
        excluded_factor_sources = c("rough_liquid_fallback"),
        drop_groups = character(0),
        group_winsor_p = 0
      )
    )

    run_robustness_scenario <- function(arguments) {
      message("Running robustness scenario: ", arguments$scenario_name)
      group_data <- tryCatch(
        do.call(
          build_group_dataset,
          c(
            list(
              item_data = item_data,
              household_data = household_data,
              group_mapping = empirical_results$group_mapping
            ),
            arguments
          )
        ),
        error = function(e) e
      )

      if (inherits(group_data, "error")) {
        warning(
          "Robustness scenario failed and was recorded rather than terminating the pipeline: ",
          arguments$scenario_name, " -- ", conditionMessage(group_data)
        )
        return(list(
          metadata = tibble(
            scenario = arguments$scenario_name,
            group_set = arguments$group_set,
            minimum_conversion_coverage = arguments$min_conversion_coverage,
            market_mode = arguments$market_mode,
            allowed_factor_sources = if (is.null(arguments$allowed_factor_sources)) {
              "all baseline matched sources"
            } else {
              paste(arguments$allowed_factor_sources, collapse = ";")
            },
            excluded_factor_sources = paste(arguments$excluded_factor_sources, collapse = ";"),
            dropped_groups = paste(arguments$drop_groups, collapse = ";"),
            group_winsor_p = arguments$group_winsor_p,
            n_households = NA_integer_,
            n_groups = NA_integer_,
            omitted_group = NA_character_,
            maximum_group_identity_error = NA_real_,
            status = "failed",
            error_message = conditionMessage(group_data)
          ),
          index = tibble(),
          expenditure_summary = tibble(),
          price_summary = tibble(),
          household_descriptives = tibble(),
          expenditure_levels = tibble(),
          price_levels = tibble(),
          fill_levels = tibble(),
          winsorization_thresholds = tibble(),
          construction_manifest = tibble(),
          household_sample_audit = tibble()
        ))
      }

      wide <- group_data$wide
      groups_safe <- group_data$groups_safe
      omitted <- groups_safe[length(groups_safe)]

      specifications <- tribble(
        ~spec,                    ~price_reg, ~price_index, ~index_type,
        "LA_PM_UVstd_current",     "market",   "uv_std",     "current",
        "LA_PM_UVraw_current",     "market",   "uv_raw",     "current"
      )

      fits <- tryCatch(
        pmap(
          specifications,
          function(spec, price_reg, price_index, index_type) {
            fit_laaids(
              data = wide,
              groups_safe = groups_safe,
              spec_name = spec,
              price_reg = price_reg,
              price_index = price_index,
              index_type = index_type,
              omitted_group = omitted,
              use_controls = TRUE
            )
          }
        ),
        error = function(e) e
      )

      if (inherits(fits, "error")) {
        warning(
          "Robustness estimation failed and was recorded rather than terminating the pipeline: ",
          arguments$scenario_name, " -- ", conditionMessage(fits),
          call. = FALSE
        )
        return(list(
          metadata = tibble(
            scenario = arguments$scenario_name,
            group_set = arguments$group_set,
            minimum_conversion_coverage = arguments$min_conversion_coverage,
            market_mode = arguments$market_mode,
            allowed_factor_sources = if (is.null(arguments$allowed_factor_sources)) {
              "all baseline matched sources"
            } else {
              paste(arguments$allowed_factor_sources, collapse = ";")
            },
            excluded_factor_sources = paste(arguments$excluded_factor_sources, collapse = ";"),
            dropped_groups = paste(arguments$drop_groups, collapse = ";"),
            group_winsor_p = arguments$group_winsor_p,
            n_households = nrow(wide),
            n_groups = length(groups_safe),
            omitted_group = omitted,
            maximum_group_identity_error = group_data$max_group_identity_error,
            status = "failed_estimation",
            error_message = conditionMessage(fits)
          ),
          index = tibble(),
          expenditure_summary = tibble(),
          price_summary = tibble(),
          household_descriptives = tibble(),
          expenditure_levels = tibble(),
          price_levels = tibble(),
          fill_levels = group_data$fill_level_summary,
          winsorization_thresholds = group_data$group_winsorization_thresholds,
          construction_manifest = group_data$construction_manifest,
          household_sample_audit = group_data$household_sample_audit
        ))
      }

      expenditure <- map_dfr(fits, "expenditure_elasticities")
      price <- map_dfr(fits, "price_elasticities")
      expenditure_cells <- compare_expenditure_specs(
        expenditure,
        "LA_PM_UVraw_current", "LA_PM_UVstd_current",
        "current_stone_only_raw_vs_standardized"
      )
      price_cells <- compare_price_specs(
        price,
        "LA_PM_UVraw_current", "LA_PM_UVstd_current",
        "current_stone_only_raw_vs_standardized"
      )

      W <- matrix_from_wide(wide, "w_", groups_safe)
      L_raw <- matrix_from_wide(wide, "lp_uv_raw_", groups_safe)
      L_std <- matrix_from_wide(wide, "lp_uv_std_", groups_safe)
      weights <- ifelse(is.finite(wide$hh_wgt) & wide$hh_wgt > 0, wide$hh_wgt, 1)

      index_results <- map_dfr(c("current", "corrected", "tornqvist", "base"), function(index_type) {
        raw_index <- make_price_index(L_raw, W, index_type, weights)$index
        std_index <- make_price_index(L_std, W, index_type, weights)$index
        difference <- raw_index - std_index
        summarise_index_difference(
          difference = difference,
          weights = weights,
          all_food_C_centered = wide$C_h_all_centered,
          system_C_centered = wide$C_h_system_centered,
          index_type = index_type,
          scenario = arguments$scenario_name
        ) %>%
          mutate(
            n_groups = length(groups_safe),
            primary_dispersion_metric =
              "mean_absolute_centered_difference"
          )
      })

      expenditure_summary <- expenditure_cells %>%
        group_by(expenditure_elasticity_definition) %>%
        summarise(
          scenario = arguments$scenario_name,
          n_households = nrow(wide),
          n_groups = length(groups_safe),
          mean_absolute_beta_difference = mean(abs(beta_raw_minus_std), na.rm = TRUE),
          maximum_absolute_beta_difference = max(abs(beta_raw_minus_std), na.rm = TRUE),
          mean_absolute_expenditure_elasticity_difference =
            mean(abs(expenditure_elasticity_raw_minus_std), na.rm = TRUE),
          maximum_absolute_expenditure_elasticity_difference =
            max(abs(expenditure_elasticity_raw_minus_std), na.rm = TRUE),
          .groups = "drop"
        )

      price_summary <- price_cells %>%
        group_by(elasticity_definition, elasticity_type) %>%
        summarise(
          scenario = arguments$scenario_name,
          n_households = nrow(wide),
          n_groups = length(groups_safe),
          mean_absolute_difference = mean(absolute_difference, na.rm = TRUE),
          median_absolute_difference = median(absolute_difference, na.rm = TRUE),
          maximum_absolute_difference = max(absolute_difference, na.rm = TRUE),
          sign_change_rate = mean(sign_change, na.rm = TRUE),
          .groups = "drop"
        )

      household_descriptives <- wide %>%
        summarise(
          scenario = arguments$scenario_name,
          n_households = n(),
          mean_conversion_coverage = mean(scenario_conversion_coverage, na.rm = TRUE),
          mean_C_h_all = mean(C_h_all, na.rm = TRUE),
          sd_C_h_all = sd(C_h_all, na.rm = TRUE),
          mean_abs_C_h_all_centered = mean(abs(C_h_all_centered), na.rm = TRUE),
          survey_weighted_mean_abs_C_h_all_centered =
            wtd_mean(abs(C_h_all - wtd_mean(C_h_all, hh_wgt)), hh_wgt),
          mean_C_h_system = mean(C_h_system, na.rm = TRUE),
          sd_C_h_system = sd(C_h_system, na.rm = TRUE),
          mean_abs_C_h_system_centered = mean(abs(C_h_system_centered), na.rm = TRUE),
          survey_weighted_mean_abs_C_h_system_centered =
            wtd_mean(abs(C_h_system - wtd_mean(C_h_system, hh_wgt)), hh_wgt),
          mean_nonkg_expenditure_share = mean(share_nonkg_exp_matched, na.rm = TRUE),
          survey_weighted_mean_nonkg_expenditure_share =
            wtd_mean(share_nonkg_exp_matched, hh_wgt)
        )

      scenario_metadata <- tibble(
        scenario = arguments$scenario_name,
        group_set = arguments$group_set,
        minimum_conversion_coverage = arguments$min_conversion_coverage,
        market_mode = arguments$market_mode,
        allowed_factor_sources = if (is.null(arguments$allowed_factor_sources)) {
          "all baseline matched sources"
        } else {
          paste(arguments$allowed_factor_sources, collapse = ";")
        },
        excluded_factor_sources = paste(arguments$excluded_factor_sources, collapse = ";"),
        dropped_groups = paste(arguments$drop_groups, collapse = ";"),
        group_winsor_p = arguments$group_winsor_p,
        n_households = nrow(wide),
        n_groups = length(groups_safe),
        omitted_group = omitted,
        maximum_group_identity_error = group_data$max_group_identity_error,
        status = "completed",
        error_message = NA_character_
      )

      scenario_dir <- file.path(config$robustness_dir, arguments$scenario_name)
      if (!dir.exists(scenario_dir)) dir.create(scenario_dir, recursive = TRUE)
      saveRDS(wide, file.path(scenario_dir, "group_wide.rds"))
      write_csv(group_data$group_summary, file.path(scenario_dir, "group_summary.csv"))
      write_csv(
        group_data$household_sample_audit,
        file.path(scenario_dir, "household_model_input_sample_audit.csv")
      )
      write_csv(group_data$market_match_summary, file.path(scenario_dir, "market_match_summary.csv"))
      write_csv(group_data$fill_level_summary, file.path(scenario_dir, "group_fill_level_summary.csv"))
      write_csv(
        group_data$group_winsorization_thresholds,
        file.path(scenario_dir, "group_winsorization_thresholds.csv")
      )
      write_csv(group_data$base_item_weights, file.path(scenario_dir, "market_group_base_item_weights.csv"))
      write_csv(
        group_data$construction_manifest,
        file.path(scenario_dir, "group_data_construction_manifest.csv")
      )
      write_csv(expenditure, file.path(scenario_dir, "expenditure_elasticity_levels.csv"))
      write_csv(price, file.path(scenario_dir, "price_elasticity_levels.csv"))
      write_csv(expenditure_cells, file.path(scenario_dir, "expenditure_raw_vs_std_cells.csv"))
      write_csv(price_cells, file.path(scenario_dir, "price_raw_vs_std_cells.csv"))

      list(
        metadata = scenario_metadata,
        index = index_results,
        expenditure_summary = expenditure_summary,
        price_summary = price_summary,
        household_descriptives = household_descriptives,
        expenditure_levels = expenditure,
        price_levels = price,
        fill_levels = group_data$fill_level_summary,
        winsorization_thresholds = group_data$group_winsorization_thresholds,
        construction_manifest = group_data$construction_manifest,
        household_sample_audit = group_data$household_sample_audit
      )
    }

    robustness_results <- purrr::map(robustness_scenarios, run_robustness_scenario)
    robustness_metadata <- map_dfr(robustness_results, "metadata")
    robustness_index <- map_dfr(robustness_results, "index")
    robustness_expenditure <- map_dfr(robustness_results, "expenditure_summary")
    robustness_price <- map_dfr(robustness_results, "price_summary")
    robustness_households <- map_dfr(robustness_results, "household_descriptives")
    robustness_fill_levels <- map_dfr(robustness_results, "fill_levels")
    robustness_winsorization_thresholds <- map_dfr(
      robustness_results, "winsorization_thresholds"
    )
    robustness_construction_manifest <- map_dfr(
      robustness_results, "construction_manifest"
    )
    robustness_household_sample_audit <- map_dfr(
      robustness_results, "household_sample_audit"
    )

    write_csv(robustness_metadata, file.path(config$robustness_dir, "robustness_scenario_metadata.csv"))
    write_csv(robustness_index, file.path(config$robustness_dir, "robustness_index_sensitivity.csv"))
    write_csv(robustness_expenditure, file.path(config$robustness_dir, "robustness_expenditure_sensitivity.csv"))
    write_csv(robustness_price, file.path(config$robustness_dir, "robustness_price_elasticity_sensitivity.csv"))
    write_csv(robustness_households, file.path(config$robustness_dir, "robustness_household_descriptives.csv"))
    write_csv(robustness_fill_levels, file.path(config$robustness_dir, "robustness_group_fill_levels.csv"))
    write_csv(
      robustness_winsorization_thresholds,
      file.path(config$robustness_dir, "robustness_group_winsorization_thresholds.csv")
    )
    write_csv(
      robustness_construction_manifest,
      file.path(config$robustness_dir, "robustness_group_construction_manifest.csv")
    )
    write_csv(
      robustness_household_sample_audit,
      file.path(config$robustness_dir, "robustness_household_model_input_sample_audit.csv")
    )

    figure_theme <- function(base_size = 11) {
      theme_bw(base_size = base_size) +
        theme(
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          strip.background = element_rect(fill = "grey90", colour = "grey60"),
          axis.title = element_text(face = "plain")
        )
    }

    save_figure <- function(plot, filename, width = 7.2, height = 5.0) {
      pdf_path <- file.path(config$figure_dir, paste0(filename, ".pdf"))
      ggsave(pdf_path, plot, width = width, height = height, units = "in")
      if (isTRUE(config$save_png)) {
        ggsave(
          file.path(config$figure_dir, paste0(filename, ".png")),
          plot, width = width, height = height, units = "in", dpi = 300
        )
      }
      invisible(pdf_path)
    }

    group_labels <- c(
      animal_protein = "Animal protein",
      cereals = "Cereals",
      oils_sugar_condiments = "Oils, sugar,\ncondiments",
      pulses_nuts = "Pulses and nuts",
      roots_tubers = "Roots and tubers",
      vegetables = "Vegetables",
      fruits = "Fruits",
      dairy = "Dairy"
    )

    baseline_wide <- empirical_results$baseline_group_data$wide

    fig2a_plot <- baseline_wide %>%
      filter(is.finite(C_h_system_centered)) %>%
      ggplot(aes(x = C_h_system_centered)) +
      geom_histogram(bins = 60, fill = "grey70", colour = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
      labs(x = expression(C[h]^system - bar(C)^system), y = "Number of households") +
      figure_theme()
    save_figure(fig2a_plot, "fig2a")

    fig2b_plot <- baseline_wide %>%
      filter(
        is.finite(share_nonkg_exp_matched),
        share_nonkg_exp_matched >= 0,
        share_nonkg_exp_matched <= 1
      ) %>%
      ggplot(aes(x = share_nonkg_exp_matched)) +
      geom_histogram(bins = 40, fill = "grey70", colour = "white") +
      geom_vline(
        xintercept = mean(baseline_wide$share_nonkg_exp_matched, na.rm = TRUE),
        linetype = "dashed", linewidth = 0.5
      ) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(
        x = "Share of conversion-matched purchased food expenditure reported in non-kg units",
        y = "Number of households"
      ) +
      figure_theme()
    save_figure(fig2b_plot, "fig2b")

    benchmark_plot <- empirical_results$benchmark_main_summary %>%
      select(
        unit_value_measure,
        expenditure_weighted_correlation,
        expenditure_weighted_mae,
        expenditure_weighted_rmse
      ) %>%
      pivot_longer(
        cols = starts_with("expenditure_weighted"),
        names_to = "statistic", values_to = "value"
      ) %>%
      mutate(
        unit_value_measure = recode(
          unit_value_measure,
          raw = "Raw unit value",
          standardized = "Standardized unit value"
        ),
        statistic = recode(
          statistic,
          expenditure_weighted_correlation = "Correlation with market-survey log price",
          expenditure_weighted_mae = "Mean absolute log-price gap",
          expenditure_weighted_rmse = "Root mean squared log-price gap"
        )
      )
    fig3_plot <- benchmark_plot %>%
      ggplot(aes(x = unit_value_measure, y = value)) +
      geom_col(fill = "grey70", colour = "grey30", width = 0.65) +
      facet_wrap(~ statistic, scales = "free_y") +
      labs(x = NULL, y = NULL) +
      figure_theme() +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
    save_figure(fig3_plot, "fig3", width = 8.0, height = 4.8)

    index_labels <- c(
      current = "Current-share Stone",
      corrected = "Corrected Stone",
      tornqvist = "T\u00f6rnqvist",
      base = "Base-share index"
    )
    fig4_point_data <- empirical_results$index_sensitivity %>%
      mutate(
        index_label = recode(index_type, !!!index_labels),
        index_label = factor(index_label, levels = rev(unname(index_labels)))
      ) %>%
      select(
        index_label,
        mean_absolute_centered_difference,
        correlation_with_all_food_C_h_centered
      ) %>%
      pivot_longer(
        cols = c(
          mean_absolute_centered_difference,
          correlation_with_all_food_C_h_centered
        ),
        names_to = "statistic", values_to = "value"
      ) %>%
      mutate(
        statistic = recode(
          statistic,
          mean_absolute_centered_difference =
            "Mean absolute centered index difference",
          correlation_with_all_food_C_h_centered =
            "Correlation with centered all-food scale component"
        )
      )

    fig4_interval_data <- tibble()
    if (!is.null(empirical_results$empirical_bootstrap_summary) &&
        !is.null(empirical_results$empirical_bootstrap_summary$index_sensitivity)) {
      fig4_interval_data <-
        empirical_results$empirical_bootstrap_summary$index_sensitivity %>%
        filter(
          metric %in% c(
            "mean_absolute_centered_difference",
            "correlation_with_all_food_C_h_centered"
          )
        ) %>%
        mutate(
          index_label = recode(index_type, !!!index_labels),
          index_label = factor(index_label, levels = rev(unname(index_labels))),
          statistic = recode(
            metric,
            mean_absolute_centered_difference =
              "Mean absolute centered index difference",
            correlation_with_all_food_C_h_centered =
              "Correlation with centered all-food scale component"
          )
        )
    }

    fig4_plot <- ggplot(
      fig4_point_data,
      aes(x = index_label, y = value)
    ) +
      geom_col(fill = "grey70", colour = "grey30", width = 0.7)
    if (nrow(fig4_interval_data) > 0) {
      fig4_plot <- fig4_plot +
        geom_errorbar(
          data = fig4_interval_data,
          aes(x = index_label, ymin = ci_lower, ymax = ci_upper),
          inherit.aes = FALSE,
          width = 0.18,
          linewidth = 0.55
        )
    }
    fig4_plot <- fig4_plot +
      coord_flip() +
      facet_wrap(~ statistic, scales = "free_y", ncol = 1) +
      labs(x = NULL, y = NULL) +
      figure_theme()
    save_figure(fig4_plot, "fig4", height = 6.2)

    expenditure_plot_data <- empirical_results$central_expenditure_levels %>%
      filter(
        expenditure_elasticity_definition == "conventional_1_plus_beta_over_share"
      )

    bootstrap_level_summary <- empirical_results$empirical_bootstrap_summary
    if (!is.null(bootstrap_level_summary)) {
      expenditure_plot_data <- expenditure_plot_data %>%
        left_join(
          bootstrap_level_summary$expenditure_levels %>%
            select(
              spec, group, expenditure_elasticity_definition,
              expenditure_elasticity_ci_lower,
              expenditure_elasticity_ci_upper
            ),
          by = c("spec", "group", "expenditure_elasticity_definition")
        )
    } else {
      expenditure_plot_data$expenditure_elasticity_ci_lower <- NA_real_
      expenditure_plot_data$expenditure_elasticity_ci_upper <- NA_real_
    }

    expenditure_plot_data <- expenditure_plot_data %>%
      mutate(
        specification = recode(
          spec,
          LA_PM_UVstd_current = "Standardized unit values in Stone index",
          LA_PM_UVraw_current = "Raw unit values in Stone index"
        ),
        group_label = recode(group, !!!group_labels, .default = group)
      )

    figure5_group_levels <- unname(group_labels[
      empirical_results$baseline_group_data$groups_safe
    ])
    figure5_group_levels <- figure5_group_levels[!is.na(figure5_group_levels)]
    expenditure_plot_data <- expenditure_plot_data %>%
      mutate(
        group_y = match(group_label, figure5_group_levels),
        specification_offset = if_else(spec == "LA_PM_UVstd_current", -0.12, 0.12),
        plot_y = group_y + specification_offset
      )

    figure5_pairs <- expenditure_plot_data %>%
      select(group, group_label, group_y, spec, expenditure_elasticity) %>%
      pivot_wider(names_from = spec, values_from = expenditure_elasticity) %>%
      filter(
        is.finite(LA_PM_UVstd_current),
        is.finite(LA_PM_UVraw_current)
      )

    fig5_plot <- ggplot() +
      geom_segment(
        data = figure5_pairs,
        aes(
          x = LA_PM_UVstd_current,
          xend = LA_PM_UVraw_current,
          y = group_y,
          yend = group_y
        ),
        colour = "grey60",
        linewidth = 0.7
      )

    if (any(
      is.finite(expenditure_plot_data$expenditure_elasticity_ci_lower) &
        is.finite(expenditure_plot_data$expenditure_elasticity_ci_upper)
    )) {
      fig5_plot <- fig5_plot +
        geom_segment(
          data = expenditure_plot_data,
          aes(
            x = expenditure_elasticity_ci_lower,
            xend = expenditure_elasticity_ci_upper,
            y = plot_y,
            yend = plot_y
          ),
          linewidth = 0.55
        )
    }

    fig5_plot <- fig5_plot +
      geom_point(
        data = expenditure_plot_data,
        aes(x = expenditure_elasticity, y = plot_y, shape = specification),
        size = 2.8
      ) +
      geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5) +
      scale_y_continuous(
        breaks = seq_along(figure5_group_levels),
        labels = figure5_group_levels
      ) +
      labs(x = "Conditional expenditure elasticity", y = NULL) +
      figure_theme()
    save_figure(fig5_plot, "fig5")

    plot_elasticity_heatmap <- function(data, filename, value_column = "elasticity") {
      group_order <- c(
        "cereals", "roots_tubers", "pulses_nuts", "vegetables",
        "animal_protein", "oils_sugar_condiments"
      )
      plot_data <- data %>%
        mutate(
          i_label = recode(i, !!!group_labels),
          j_label = recode(j, !!!group_labels),
          i_label = factor(i_label, levels = recode(group_order, !!!group_labels)),
          j_label = factor(j_label, levels = recode(group_order, !!!group_labels))
        )
      limit <- max(abs(plot_data[[value_column]]), na.rm = TRUE)
      if (!is.finite(limit) || limit == 0) limit <- 1
      plot <- ggplot(plot_data, aes(x = j_label, y = forcats::fct_rev(i_label), fill = .data[[value_column]])) +
        geom_tile(colour = "white") +
        geom_text(aes(label = sprintf("%.2f", .data[[value_column]])), size = 3) +
        scale_fill_gradient2(
          low = "grey15", mid = "white", high = "grey75",
          midpoint = 0, limits = c(-limit, limit)
        ) +
        labs(x = "Price group j", y = "Demand group i", fill = NULL) +
        figure_theme() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      save_figure(plot, filename, width = 7.2, height = 5.8)
    }

    supplement_figure_files <- list(
      partial_regressor_index_fixed = c(
        standardized = "figS1",
        raw = "figS2",
        difference = "figS3"
      ),
      joint_one_for_one_Green_Alston = c(
        standardized = "figS4",
        raw = "figS5",
        difference = "figS6"
      )
    )

    for (definition in names(supplement_figure_files)) {
      filenames <- supplement_figure_files[[definition]]
      standardized <- empirical_results$central_price_levels %>%
        filter(spec == "LA_PM_UVstd_current", elasticity_definition == definition)
      raw <- empirical_results$central_price_levels %>%
        filter(spec == "LA_PM_UVraw_current", elasticity_definition == definition)
      difference <- raw %>%
        select(i, j, elasticity_raw = elasticity) %>%
        inner_join(
          standardized %>% select(i, j, elasticity_standardized = elasticity),
          by = c("i", "j")
        ) %>%
        mutate(elasticity_difference = elasticity_raw - elasticity_standardized)

      plot_elasticity_heatmap(
        standardized,
        unname(filenames[["standardized"]])
      )
      plot_elasticity_heatmap(
        raw,
        unname(filenames[["raw"]])
      )
      plot_elasticity_heatmap(
        difference,
        unname(filenames[["difference"]]),
        value_column = "elasticity_difference"
      )
    }

    if (any(robustness_metadata$status != "completed")) {
      stop("One or more robustness scenarios failed; inspect robustness_scenario_metadata.csv.", call. = FALSE)
    }
    invisible(robustness_metadata)
  }
  message("Empirical sample: ", nrow(wide_baseline), " households. Generating robustness results and figures.")
  robustness_metadata <- run_robustness_and_figures()
  invisible(list(sample_size = nrow(wide_baseline),
                 bootstrap_design = empirical_bootstrap$design,
                 robustness_metadata = robustness_metadata))
}

run_monte_carlo <- function(config) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for Monte Carlo figures.", call. = FALSE)
  }
  output_dir <- config$mc_output_dir
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  mc_raw_results_path <- file.path(output_dir, "mc_raw_results.rds")
  mc_configuration_provenance <- "executed_configuration_embedded_in_raw_rds"
  mc_true_parameter_provenance <- "generated_from_executed_configuration"
  mc_rng_kind_recorded <- paste(RNGkind(), collapse = ";")

  # Exact AIDS and price indices
  safe_rbind <- function(x) {
    x <- x[!vapply(x, is.null, logical(1))]
    if (length(x) == 0) return(NULL)
    do.call(rbind, x)
  }

  mc_numerical_tol <- config$mc_numerical_tol
  mc_optimizer_maxit <- 1500L
  mc_optimizer_reltol <- 1e-10
  mc_retry_shift_scale <- 0.01

  share_above_tolerance <- function(x, tol = mc_numerical_tol) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    mean(x > tol)
  }

  sign_with_tolerance <- function(x, tol = mc_numerical_tol) {
    out <- sign(x)
    out[is.finite(x) & abs(x) <= tol] <- 0
    out
  }

  check_price_corr <- function(rho, n) {
    lower <- -1 / (n - 1)
    if (rho <= lower || rho >= 1) {
      stop(
        "For an equicorrelation matrix with n goods, price_corr must be in (",
        round(lower, 4), ", 1)."
      )
    }
  }

  complete_gamma_from_block <- function(H) {
    k <- nrow(H)
    n <- k + 1

    G <- matrix(0, n, n)
    G[seq_len(k), seq_len(k)] <- H

    for (i in seq_len(k)) {
      G[i, n] <- -sum(G[i, seq_len(k)])
    }
    for (j in seq_len(k)) {
      G[n, j] <- -sum(G[seq_len(k), j])
    }
    G[n, n] <- -sum(G[n, seq_len(k)])

    G
  }

  pack_theta <- function(par) {
    n <- length(par$alpha)
    k <- n - 1
    H <- par$Gamma[seq_len(k), seq_len(k), drop = FALSE]
    idx <- upper.tri(H, diag = TRUE)

    c(
      par$alpha[seq_len(k)],
      par$beta[seq_len(k)],
      H[idx]
    )
  }

  unpack_theta <- function(theta, n, alpha0 = 0) {
    k <- n - 1

    alpha_free <- theta[seq_len(k)]
    beta_free  <- theta[k + seq_len(k)]

    alpha <- c(alpha_free, 1 - sum(alpha_free))
    beta  <- c(beta_free, -sum(beta_free))

    n_H <- k * (k + 1) / 2
    H_vals <- theta[2 * k + seq_len(n_H)]

    H <- matrix(0, k, k)
    idx <- upper.tri(H, diag = TRUE)
    H[idx] <- H_vals
    H[lower.tri(H)] <- t(H)[lower.tri(H)]

    Gamma <- complete_gamma_from_block(H)

    list(
      n = n,
      alpha0 = alpha0,
      alpha = alpha,
      beta = beta,
      Gamma = Gamma
    )
  }

  aids_index <- function(L, par) {
    as.vector(
      par$alpha0 +
        L %*% par$alpha +
        0.5 * rowSums((L %*% par$Gamma) * L)
    )
  }

  aids_price_index_weights <- function(L, par) {
    matrix(par$alpha, nrow = nrow(L), ncol = ncol(L), byrow = TRUE) +
      L %*% t(par$Gamma)
  }

  transform_aids_units <- function(par, c_vec) {
    par_star <- par
    par_star$alpha0 <-
      par$alpha0 -
      sum(par$alpha * c_vec) +
      0.5 * as.numeric(t(c_vec) %*% par$Gamma %*% c_vec)
    par_star$alpha <- par$alpha - as.vector(par$Gamma %*% c_vec)
    par_star
  }

  predict_shares <- function(par, L, m, exact = TRUE, index = NULL) {
    N <- nrow(L)
    n <- ncol(L)

    if (exact) {
      real_exp <- m - aids_index(L, par)
    } else {
      if (is.null(index)) {
        stop("For LA/AIDS, an approximating price index must be supplied.")
      }
      real_exp <- m - index
    }

    matrix(par$alpha, nrow = N, ncol = n, byrow = TRUE) +
      L %*% t(par$Gamma) +
      real_exp %o% par$beta
  }

  make_true_params <- function(beta = c(0.04, -0.02, 0.03, -0.05),
                               gamma_scale = 1,
                               homothetic = FALSE,
                               alpha0 = 0) {

    alpha <- c(0.35, 0.25, 0.25, 0.15)

    if (homothetic) {
      beta <- rep(0, 4)
    } else {
      if (abs(sum(beta)) > 1e-10) {
        stop("The supplied beta vector must sum to zero.")
      }
    }

    H_base <- matrix(
      c(
        -0.020,  0.006,  0.004,
         0.006, -0.018,  0.005,
         0.004,  0.005, -0.015
      ),
      nrow = 3,
      byrow = TRUE
    )

    H <- gamma_scale * H_base
    Gamma <- complete_gamma_from_block(H)

    list(
      n = 4,
      alpha0 = alpha0,
      alpha = alpha,
      beta = beta,
      Gamma = Gamma
    )
  }

  make_price_index <- function(L, W,
                               type = "stone_current",
                               base_w = NULL,
                               base_l = NULL) {

    N <- nrow(L)
    n <- ncol(L)

    if (is.null(base_w)) base_w <- colMeans(W)
    if (is.null(base_l)) base_l <- colMeans(L)

    base_w_mat <- matrix(base_w, nrow = N, ncol = n, byrow = TRUE)
    base_l_mat <- matrix(base_l, nrow = N, ncol = n, byrow = TRUE)

    if (type == "stone_current") {
      return(rowSums(W * L))
    }

    if (type == "stone_base") {
      return(as.vector(L %*% base_w))
    }

    if (type == "laspeyres") {
      return(as.vector((L - base_l_mat) %*% base_w))
    }

    if (type == "paasche" || type == "corrected_stone") {
      return(rowSums(W * (L - base_l_mat)))
    }

    if (type == "tornqvist") {
      return(rowSums(0.5 * (W + base_w_mat) * (L - base_l_mat)))
    }

    stop("Unknown price index type.")
  }

  index_derivative_components <- function(index_type,
                                          L_index_eval,
                                          W_eval,
                                          base_w = NULL,
                                          base_l = NULL) {
    n <- length(W_eval)
    if (is.null(base_w)) base_w <- W_eval
    if (is.null(base_l)) base_l <- rep(0, n)

    if (index_type == "stone_current") {
      direct <- W_eval
      cvec <- W_eval * L_index_eval
      return(list(direct = direct, cvec = cvec))
    }

    if (index_type == "paasche" || index_type == "corrected_stone") {
      Ldiff <- L_index_eval - base_l
      direct <- W_eval
      cvec <- W_eval * Ldiff
      return(list(direct = direct, cvec = cvec))
    }

    if (index_type == "stone_base" || index_type == "laspeyres") {
      direct <- base_w
      cvec <- rep(0, n)
      return(list(direct = direct, cvec = cvec))
    }

    if (index_type == "tornqvist") {
      Ldiff <- L_index_eval - base_l
      direct <- 0.5 * (W_eval + base_w)
      cvec <- 0.5 * W_eval * Ldiff
      return(list(direct = direct, cvec = cvec))
    }

    stop("Unknown index_type in derivative components.")
  }

  # Data generation
  generate_data <- function(par, cfg) {
    n <- par$n
    N <- cfg$N
    check_price_corr(cfg$price_corr, n)

    for (attempt in seq_len(cfg$max_tries)) {

      Sigma_p <- cfg$sigma_logp^2 *
        ((1 - cfg$price_corr) * diag(n) +
           cfg$price_corr * matrix(1, n, n))

      L_true <- matrix(rnorm(N * n), N, n) %*% chol(Sigma_p)

      M <- rnorm(N, mean = cfg$M_mean, sd = cfg$M_sd)

      U_raw <- matrix(rnorm(N * n, sd = cfg$sigma_u), N, n)
      U <- U_raw - rowMeans(U_raw)

      W_true <-
        matrix(par$alpha, nrow = N, ncol = n, byrow = TRUE) +
        L_true %*% t(par$Gamma) +
        M %o% par$beta +
        U

      if (min(W_true) > cfg$min_share && max(W_true) < 1 - cfg$min_share) {
        break
      }

      if (attempt == cfg$max_tries) {
        stop("Could not generate valid shares. Try smaller variances.")
      }
    }

    a_true <- aids_index(L_true, par)
    m_true <- M + a_true

    quality_income <- cfg$quality_income
    quality_price  <- cfg$quality_price
    taste_loading  <- cfg$taste_loading

    Tau_income <- M %o% quality_income

    Tau_price <- -sweep(L_true, 2, quality_price, `*`)

    Tau_taste <- taste_loading * U

    A_exp <- matrix(rnorm(N * n, sd = cfg$sigma_exp_error), N, n)
    B_qty <- matrix(rnorm(N * n, sd = cfg$sigma_qty_error), N, n)
    Tau_reporting <- A_exp - B_qty

    n_clusters <- cfg$n_clusters
    cluster_id <- sample(seq_len(n_clusters), N, replace = TRUE)
    K_cluster_draws <- matrix(
      rnorm(n_clusters * n, sd = cfg$sigma_local_cluster),
      n_clusters,
      n
    )
    K_cluster <- K_cluster_draws[cluster_id, , drop = FALSE]
    K_household <- matrix(rnorm(N * n, sd = cfg$sigma_local_household), N, n)
    Kappa <- K_cluster + K_household

    Tau_extra <- matrix(rnorm(N * n, sd = cfg$sigma_extra_uv), N, n)

    Tau <- Tau_income + Tau_price + Tau_taste + Tau_reporting + Kappa + Tau_extra

    c_vec <- cfg$unit_scale
    if (length(c_vec) != n) stop("unit_scale must have length equal to number of goods.")
    C_mat <- matrix(c_vec, nrow = N, ncol = n, byrow = TRUE)

    L_uv_norm <- L_true + Tau          # fixed unit scale removed
    L_uv_raw  <- L_true + C_mat + Tau  # raw household unit values

    if (isTRUE(cfg$share_reporting_error)) {
      A_for_shares <- A_exp
    } else {
      A_for_shares <- matrix(0, N, n)
    }

    exp_factor <- rowSums(W_true * exp(A_for_shares))
    W_obs <- W_true * exp(A_for_shares) / exp_factor
    m_obs <- m_true + log(exp_factor)

    C_h_obs <- as.vector(W_obs %*% c_vec)
    C_h_true <- as.vector(W_true %*% c_vec)

    list(
      L_true = L_true,
      L_uv_norm = L_uv_norm,
      L_uv_raw = L_uv_raw,
      W_true = W_true,
      W_obs = W_obs,
      m_true = m_true,
      m_obs = m_obs,
      M = M,
      U = U,
      Tau = Tau,
      Tau_income = Tau_income,
      Tau_price = Tau_price,
      Tau_taste = Tau_taste,
      Tau_reporting = Tau_reporting,
      Kappa = Kappa,
      A_exp = A_exp,
      B_qty = B_qty,
      unit_scale = c_vec,
      C_h_obs = C_h_obs,
      C_h_true = C_h_true,
      cluster_id = cluster_id
    )
  }

  # Estimation
  make_start_from_la <- function(W, L, m, index = NULL) {
    n <- ncol(W)
    k <- n - 1

    if (is.null(index)) {
      index <- rowSums(W * L)
    }

    real_exp <- m - index
    X <- cbind(1, L, real_exp)

    alpha_start <- rep(NA_real_, k)
    beta_start  <- rep(NA_real_, k)
    gamma_rows  <- matrix(NA_real_, k, n)

    for (i in seq_len(k)) {
      fit_i <- lm.fit(X, W[, i])
      coef_i <- fit_i$coefficients
      coef_i[is.na(coef_i)] <- 0

      alpha_start[i] <- coef_i[1]
      gamma_rows[i, ] <- coef_i[2:(n + 1)]
      beta_start[i] <- coef_i[n + 2]
    }

    H <- gamma_rows[, seq_len(k), drop = FALSE]
    H <- 0.5 * (H + t(H))

    Gamma_start <- complete_gamma_from_block(H)

    par_start <- list(
      n = n,
      alpha0 = 0,
      alpha = c(alpha_start, 1 - sum(alpha_start)),
      beta = c(beta_start, -sum(beta_start)),
      Gamma = Gamma_start
    )

    pack_theta(par_start)
  }

  fit_aids_model <- function(W, L, m,
                             exact = TRUE,
                             index = NULL,
                             start_theta = NULL,
                             control = list(
                               maxit = mc_optimizer_maxit,
                               reltol = mc_optimizer_reltol
                             )) {
    n <- ncol(W)
    k <- n - 1

    if (is.null(start_theta)) {
      start_index <- if (is.null(index)) rowSums(W * L) else index
      start_theta <- make_start_from_la(W, L, m, index = start_index)
    }

    obj <- function(theta) {
      par_hat <- unpack_theta(theta, n = n, alpha0 = 0)
      W_hat <- predict_shares(
        par = par_hat,
        L = L,
        m = m,
        exact = exact,
        index = index
      )
      resid <- W[, seq_len(k), drop = FALSE] - W_hat[, seq_len(k), drop = FALSE]
      value <- sum(resid^2)
      if (!is.finite(value)) .Machine$double.xmax else value
    }

    # Deterministic retries preserve the simulation random-number stream.
    deterministic_shift <- mc_retry_shift_scale * sin(seq_along(start_theta))
    attempts <- list(
      list(method = "BFGS", start = start_theta),
      list(method = "BFGS", start = start_theta + deterministic_shift),
      list(method = "BFGS", start = start_theta - deterministic_shift),
      list(method = "Nelder-Mead", start = start_theta)
    )

    best_nonconverged <- NULL
    n_attempted <- 0L
    for (a in seq_along(attempts)) {
      n_attempted <- n_attempted + 1L
      candidate <- try(
        optim(
          par = attempts[[a]]$start,
          fn = obj,
          method = attempts[[a]]$method,
          control = control
        ),
        silent = TRUE
      )
      if (inherits(candidate, "try-error")) next

      if (identical(candidate$convergence, 0L)) {
        return(list(
          ok = TRUE,
          par = unpack_theta(candidate$par, n = n, alpha0 = 0),
          value = candidate$value,
          convergence = candidate$convergence,
          message = paste0("OK after ", n_attempted, " attempt(s)")
        ))
      }

      if (is.null(best_nonconverged) || candidate$value < best_nonconverged$value) {
        best_nonconverged <- candidate
      }
    }

    if (is.null(best_nonconverged)) {
      return(list(
        ok = FALSE,
        par = NULL,
        value = NA_real_,
        convergence = NA_integer_,
        message = "all optimization attempts failed"
      ))
    }

    list(
      ok = FALSE,
      par = NULL,
      value = best_nonconverged$value,
      convergence = best_nonconverged$convergence,
      message = "no converged optimization attempt; estimate excluded"
    )
  }

  # Elasticities
  elasticity_aids_price <- function(par, L_eval, W_eval) {
    n <- length(W_eval)
    s_eval <- as.vector(par$alpha + par$Gamma %*% L_eval)
    E <- matrix(NA_real_, n, n)

    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        delta_ij <- as.numeric(i == j)
        E[i, j] <-
          -delta_ij +
          par$Gamma[i, j] / W_eval[i] -
          (par$beta[i] / W_eval[i]) * s_eval[j]
      }
    }

    E
  }

  elasticity_la_price_prime_general <- function(par, W_eval, direct) {
    n <- length(W_eval)
    E <- matrix(NA_real_, n, n)

    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        delta_ij <- as.numeric(i == j)
        E[i, j] <-
          -delta_ij +
          par$Gamma[i, j] / W_eval[i] -
          (par$beta[i] / W_eval[i]) * direct[j]
      }
    }

    E
  }

  elasticity_la_price_GA_general <- function(par, W_eval, direct, cvec) {
    n <- length(W_eval)
    E <- matrix(NA_real_, n, n)

    bvec <- par$beta / W_eval
    # Solve the joint share-weight response in Green-Alston elasticities.
    A_mat <- diag(n) + bvec %o% cvec

    for (j in seq_len(n)) {
      a_col <- rep(NA_real_, n)
      for (i in seq_len(n)) {
        delta_ij <- as.numeric(i == j)
        a_col[i] <-
          -delta_ij +
          par$Gamma[i, j] / W_eval[i] -
          bvec[i] * (direct[j] + cvec[j])
      }

      sol <- try(solve(A_mat, a_col), silent = TRUE)
      if (inherits(sol, "try-error")) {
        E[, j] <- NA_real_
      } else {
        E[, j] <- as.vector(sol)
      }
    }

    E
  }

  elasticity_aids_expenditure <- function(par, W_eval) {
    1 + par$beta / W_eval
  }

  # Specifications and records
  make_specifications <- function(dat, cfg) {
    specs <- list(
      list(
        spec = "AIDS_P",
        model = "AIDS",
        exact = TRUE,
        L_reg = dat$L_true,
        L_index = NULL,
        index_type = NA_character_,
        price_regressor = "P",
        index_price = "translog"
      ),
      list(
        spec = "AIDS_UV_norm",
        model = "AIDS",
        exact = TRUE,
        L_reg = dat$L_uv_norm,
        L_index = NULL,
        index_type = NA_character_,
        price_regressor = "UV_norm",
        index_price = "translog"
      ),
      list(
        spec = "LA_PP_current",
        model = "LA_AIDS",
        exact = FALSE,
        L_reg = dat$L_true,
        L_index = dat$L_true,
        index_type = "stone_current",
        price_regressor = "P",
        index_price = "P"
      ),
      list(
        spec = "LA_UVnorm_UVnorm_current",
        model = "LA_AIDS",
        exact = FALSE,
        L_reg = dat$L_uv_norm,
        L_index = dat$L_uv_norm,
        index_type = "stone_current",
        price_regressor = "UV_norm",
        index_price = "UV_norm"
      ),
      list(
        spec = "LA_UVraw_UVraw_current",
        model = "LA_AIDS",
        exact = FALSE,
        L_reg = dat$L_uv_raw,
        L_index = dat$L_uv_raw,
        index_type = "stone_current",
        price_regressor = "UV_raw",
        index_price = "UV_raw"
      ),
      list(
        spec = "LA_P_UVnorm_current",
        model = "LA_AIDS",
        exact = FALSE,
        L_reg = dat$L_true,
        L_index = dat$L_uv_norm,
        index_type = "stone_current",
        price_regressor = "P",
        index_price = "UV_norm"
      ),
      list(
        spec = "LA_P_UVraw_current",
        model = "LA_AIDS",
        exact = FALSE,
        L_reg = dat$L_true,
        L_index = dat$L_uv_raw,
        index_type = "stone_current",
        price_regressor = "P",
        index_price = "UV_raw"
      ),
      list(
        spec = "LA_UVnorm_P_current",
        model = "LA_AIDS",
        exact = FALSE,
        L_reg = dat$L_uv_norm,
        L_index = dat$L_true,
        index_type = "stone_current",
        price_regressor = "UV_norm",
        index_price = "P"
      ),
      list(
        spec = "LA_UVraw_P_current",
        model = "LA_AIDS",
        exact = FALSE,
        L_reg = dat$L_uv_raw,
        L_index = dat$L_true,
        index_type = "stone_current",
        price_regressor = "UV_raw",
        index_price = "P"
      )
    )

    if (isTRUE(cfg$run_invariant_indices)) {
      specs <- c(
        specs,
        list(
          list(
            spec = "LA_UVraw_UVraw_corrected",
            model = "LA_AIDS",
            exact = FALSE,
            L_reg = dat$L_uv_raw,
            L_index = dat$L_uv_raw,
            index_type = "corrected_stone",
            price_regressor = "UV_raw",
            index_price = "UV_raw"
          ),
          list(
            spec = "LA_UVraw_UVraw_tornqvist",
            model = "LA_AIDS",
            exact = FALSE,
            L_reg = dat$L_uv_raw,
            L_index = dat$L_uv_raw,
            index_type = "tornqvist",
            price_regressor = "UV_raw",
            index_price = "UV_raw"
          ),
          list(
            spec = "LA_UVnorm_UVnorm_corrected",
            model = "LA_AIDS",
            exact = FALSE,
            L_reg = dat$L_uv_norm,
            L_index = dat$L_uv_norm,
            index_type = "corrected_stone",
            price_regressor = "UV_norm",
            index_price = "UV_norm"
          ),
          list(
            spec = "LA_UVnorm_UVnorm_tornqvist",
            model = "LA_AIDS",
            exact = FALSE,
            L_reg = dat$L_uv_norm,
            L_index = dat$L_uv_norm,
            index_type = "tornqvist",
            price_regressor = "UV_norm",
            index_price = "UV_norm"
          )
        )
      )
    }

    specs
  }

  make_elasticity_records <- function(rep_id,
                                      scenario,
                                      spec_meta,
                                      formula,
                                      E_hat,
                                      E_true) {
    n <- nrow(E_true)
    grid <- expand.grid(i = seq_len(n), j = seq_len(n))

    data.frame(
      scenario = scenario,
      rep = rep_id,
      spec = spec_meta$spec,
      model = spec_meta$model,
      price_regressor = spec_meta$price_regressor,
      index_price = spec_meta$index_price,
      index_type = ifelse(is.na(spec_meta$index_type), "translog", spec_meta$index_type),
      formula = formula,
      i = grid$i,
      j = grid$j,
      type = ifelse(grid$i == grid$j, "own", "cross"),
      estimate = as.vector(E_hat),
      true = as.vector(E_true)
    )
  }

  make_beta_records <- function(rep_id, scenario, spec_meta, fit, par_true) {
    n <- par_true$n
    if (!isTRUE(fit$ok) || is.null(fit$par)) {
      beta_hat <- rep(NA_real_, n)
    } else {
      beta_hat <- fit$par$beta
    }

    data.frame(
      scenario = scenario,
      rep = rep_id,
      spec = spec_meta$spec,
      model = spec_meta$model,
      price_regressor = spec_meta$price_regressor,
      index_price = spec_meta$index_price,
      index_type = ifelse(is.na(spec_meta$index_type), "translog", spec_meta$index_type),
      i = seq_len(n),
      beta_hat = beta_hat,
      beta_true = par_true$beta,
      beta_error = beta_hat - par_true$beta
    )
  }

  make_theoretical_records <- function(rep_id, scenario, dat, par_true) {
    n <- par_true$n
    s_true <- aids_price_index_weights(dat$L_true, par_true)
    diff_ws <- dat$W_obs - s_true

    C_centered <- dat$C_h_obs - mean(dat$C_h_obs)
    K_index <- rowSums(dat$W_obs * dat$Kappa)
    tau_inner <- rowSums(diff_ws * dat$Tau)

    rows <- vector("list", n + 1L)
    for (i in seq_len(n)) {
      D_comm <- -par_true$beta[i] * C_centered
      D_tau <- -par_true$beta[i] * tau_inner
      D_total <- D_comm + D_tau
      rows[[i]] <- data.frame(
        scenario = scenario,
        rep = rep_id,
        i = i,
        mean_abs_C_centered = mean(abs(C_centered)),
        sd_C_centered = sd(C_centered),
        mean_abs_K_index = mean(abs(K_index)),
        mean_abs_tau_inner = mean(abs(tau_inner)),
        mean_abs_D_comm = mean(abs(D_comm)),
        mean_abs_D_tau = mean(abs(D_tau)),
        mean_abs_D_total = mean(abs(D_total)),
        rmse_D_comm = sqrt(mean(D_comm^2)),
        rmse_D_tau = sqrt(mean(D_tau^2)),
        rmse_D_total = sqrt(mean(D_total^2))
      )
    }

    by_good <- safe_rbind(rows[seq_len(n)])
    rows[[n + 1L]] <- data.frame(
      scenario = scenario,
      rep = rep_id,
      i = 0,
      mean_abs_C_centered = mean(by_good$mean_abs_C_centered),
      sd_C_centered = mean(by_good$sd_C_centered),
      mean_abs_K_index = mean(by_good$mean_abs_K_index),
      mean_abs_tau_inner = mean(by_good$mean_abs_tau_inner),
      mean_abs_D_comm = mean(by_good$mean_abs_D_comm),
      mean_abs_D_tau = mean(by_good$mean_abs_D_tau),
      mean_abs_D_total = mean(by_good$mean_abs_D_total),
      rmse_D_comm = mean(by_good$rmse_D_comm),
      rmse_D_tau = mean(by_good$rmse_D_tau),
      rmse_D_total = mean(by_good$rmse_D_total)
    )

    safe_rbind(rows)
  }

  # Replication runner
  run_one_replication <- function(rep_id, scenario, par_true, cfg) {
    dat <- generate_data(par_true, cfg)

    W <- dat$W_obs
    m <- dat$m_obs
    W_eval <- colMeans(W)

    W_true_eval <- colMeans(dat$W_true)
    L_true_eval <- colMeans(dat$L_true)

    E_true <- elasticity_aids_price(
      par = par_true,
      L_eval = L_true_eval,
      W_eval = W_true_eval
    )

    specs <- make_specifications(dat, cfg)

    elasticity_list <- list()
    beta_list <- list()
    convergence_list <- list()

    for (s in seq_along(specs)) {
      spec <- specs[[s]]

      if (spec$exact) {
        fit <- fit_aids_model(
          W = W,
          L = spec$L_reg,
          m = m,
          exact = TRUE
        )
      } else {
        base_w <- colMeans(W)
        base_l <- colMeans(spec$L_index)
        idx <- make_price_index(
          L = spec$L_index,
          W = W,
          type = spec$index_type,
          base_w = base_w,
          base_l = base_l
        )
        fit <- fit_aids_model(
          W = W,
          L = spec$L_reg,
          m = m,
          exact = FALSE,
          index = idx
        )
      }

      convergence_list[[length(convergence_list) + 1]] <- data.frame(
        scenario = scenario,
        rep = rep_id,
        spec = spec$spec,
        model = spec$model,
        price_regressor = spec$price_regressor,
        index_price = spec$index_price,
        index_type = ifelse(is.na(spec$index_type), "translog", spec$index_type),
        ok = fit$ok,
        convergence = fit$convergence,
        objective = fit$value,
        message = fit$message
      )

      beta_list[[length(beta_list) + 1]] <- make_beta_records(
        rep_id = rep_id,
        scenario = scenario,
        spec_meta = spec,
        fit = fit,
        par_true = par_true
      )

      if (!isTRUE(fit$ok) || is.null(fit$par)) {
        E_hat <- matrix(NA_real_, par_true$n, par_true$n)
        elasticity_list[[length(elasticity_list) + 1]] <- make_elasticity_records(
          rep_id = rep_id,
          scenario = scenario,
          spec_meta = spec,
          formula = ifelse(spec$exact, "AIDS", "Green_Alston"),
          E_hat = E_hat,
          E_true = E_true
        )
        if (!spec$exact) {
          elasticity_list[[length(elasticity_list) + 1]] <- make_elasticity_records(
            rep_id = rep_id,
            scenario = scenario,
            spec_meta = spec,
            formula = "LA_prime",
            E_hat = E_hat,
            E_true = E_true
          )
        }
        next
      }

      if (spec$exact) {
        E_hat <- elasticity_aids_price(
          par = fit$par,
          L_eval = colMeans(spec$L_reg),
          W_eval = W_eval
        )
        elasticity_list[[length(elasticity_list) + 1]] <- make_elasticity_records(
          rep_id = rep_id,
          scenario = scenario,
          spec_meta = spec,
          formula = "AIDS",
          E_hat = E_hat,
          E_true = E_true
        )
      } else {
        base_w <- colMeans(W)
        base_l <- colMeans(spec$L_index)
        comp <- index_derivative_components(
          index_type = spec$index_type,
          L_index_eval = colMeans(spec$L_index),
          W_eval = W_eval,
          base_w = base_w,
          base_l = base_l
        )

        E_hat_GA <- elasticity_la_price_GA_general(
          par = fit$par,
          W_eval = W_eval,
          direct = comp$direct,
          cvec = comp$cvec
        )
        E_hat_prime <- elasticity_la_price_prime_general(
          par = fit$par,
          W_eval = W_eval,
          direct = comp$direct
        )

        elasticity_list[[length(elasticity_list) + 1]] <- make_elasticity_records(
          rep_id = rep_id,
          scenario = scenario,
          spec_meta = spec,
          formula = "Green_Alston",
          E_hat = E_hat_GA,
          E_true = E_true
        )
        elasticity_list[[length(elasticity_list) + 1]] <- make_elasticity_records(
          rep_id = rep_id,
          scenario = scenario,
          spec_meta = spec,
          formula = "LA_prime",
          E_hat = E_hat_prime,
          E_true = E_true
        )
      }
    }

    theoretical <- make_theoretical_records(
      rep_id = rep_id,
      scenario = scenario,
      dat = dat,
      par_true = par_true
    )

    par_star <- transform_aids_units(par_true, dat$unit_scale)
    L_shifted <- sweep(dat$L_true, 2, dat$unit_scale, "+")
    max_index_gap <- max(abs(aids_index(dat$L_true, par_true) - aids_index(L_shifted, par_star)))
    oracle <- data.frame(
      scenario = scenario,
      rep = rep_id,
      max_exact_AIDS_index_gap_after_oracle_unit_transform = max_index_gap
    )

    list(
      elasticities = safe_rbind(elasticity_list),
      beta = safe_rbind(beta_list),
      convergence = safe_rbind(convergence_list),
      theoretical = theoretical,
      oracle = oracle
    )
  }

  run_monte_carlo_scenario <- function(scenario_name, cfg) {
    if (!is.null(cfg$seed)) set.seed(cfg$seed)

    par_true <- make_true_params(
      beta = cfg$beta_true,
      gamma_scale = cfg$gamma_scale,
      homothetic = cfg$homothetic,
      alpha0 = 0
    )

    elasticities_list <- vector("list", cfg$R)
    beta_list <- vector("list", cfg$R)
    convergence_list <- vector("list", cfg$R)
    theoretical_list <- vector("list", cfg$R)
    oracle_list <- vector("list", cfg$R)
    execution_list <- vector("list", cfg$R)

    for (r in seq_len(cfg$R)) {
      if (isTRUE(cfg$verbose) && (r %% cfg$print_every == 0)) {
        cat("Scenario", scenario_name, ": replication", r, "of", cfg$R, "\n")
      }

      out_r <- tryCatch(
        run_one_replication(
          rep_id = r,
          scenario = scenario_name,
          par_true = par_true,
          cfg = cfg
        ),
        error = function(e) e
      )

      if (inherits(out_r, "error")) {
        execution_list[[r]] <- data.frame(
          scenario = scenario_name,
          rep = r,
          scenario_seed = if (is.null(cfg$seed)) NA_integer_ else as.integer(cfg$seed),
          rng_stream = "sequential_within_scenario",
          status = "failed_before_model_records",
          model_specifications_attempted = 0L,
          model_specifications_converged = 0L,
          error_message = conditionMessage(out_r),
          stringsAsFactors = FALSE
        )
        warning(
          "Scenario ", scenario_name, ", replication ", r,
          " failed completely: ", conditionMessage(out_r)
        )
        next
      }

      elasticities_list[[r]] <- out_r$elasticities
      beta_list[[r]] <- out_r$beta
      convergence_list[[r]] <- out_r$convergence
      theoretical_list[[r]] <- out_r$theoretical
      oracle_list[[r]] <- out_r$oracle

      n_attempted <- if (is.null(out_r$convergence)) 0L else nrow(out_r$convergence)
      n_converged <- if (is.null(out_r$convergence)) {
        0L
      } else {
        sum(out_r$convergence$ok %in% TRUE, na.rm = TRUE)
      }
      execution_list[[r]] <- data.frame(
        scenario = scenario_name,
        rep = r,
        scenario_seed = if (is.null(cfg$seed)) NA_integer_ else as.integer(cfg$seed),
        rng_stream = "sequential_within_scenario",
        status = if (n_attempted > 0L && n_converged == n_attempted) {
          "completed_all_models_converged"
        } else {
          "completed_with_model_nonconvergence"
        },
        model_specifications_attempted = as.integer(n_attempted),
        model_specifications_converged = as.integer(n_converged),
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    list(
      scenario = scenario_name,
      configuration = cfg,
      true_parameters = par_true,
      elasticities = safe_rbind(elasticities_list),
      beta = safe_rbind(beta_list),
      convergence = safe_rbind(convergence_list),
      theoretical = safe_rbind(theoretical_list),
      oracle = safe_rbind(oracle_list),
      replication_execution = safe_rbind(execution_list)
    )
  }

  # Summary functions
  summarise_elasticities <- function(res) {
    res <- res[is.finite(res$estimate) & is.finite(res$true), , drop = FALSE]
    if (nrow(res) == 0) return(NULL)
    res$error <- res$estimate - res$true
    res$sq_error <- res$error^2
    res$abs_error <- abs(res$error)
    res$wrong_sign <- as.integer(
      sign_with_tolerance(res$estimate) != sign_with_tolerance(res$true)
    )

    key <- interaction(
      res$scenario,
      res$spec,
      res$model,
      res$price_regressor,
      res$index_price,
      res$index_type,
      res$formula,
      res$type,
      drop = TRUE
    )

    out <- lapply(split(res, key), function(z) {
      data.frame(
        scenario = z$scenario[1],
        spec = z$spec[1],
        model = z$model[1],
        price_regressor = z$price_regressor[1],
        index_price = z$index_price[1],
        index_type = z$index_type[1],
        formula = z$formula[1],
        type = z$type[1],
        n_obs = nrow(z),
        mean_true = mean(z$true),
        mean_estimate = mean(z$estimate),
        bias = mean(z$error),
        mean_abs_error = mean(z$abs_error),
        rmse = sqrt(mean(z$sq_error)),
        wrong_sign_rate = mean(z$wrong_sign)
      )
    })

    safe_rbind(out)
  }

  get_spec_estimates <- function(res, spec_name, formula_name) {
    subset(
      res,
      spec == spec_name & formula == formula_name,
      select = c("scenario", "rep", "i", "j", "type", "true", "estimate")
    )
  }

  compare_two_specs <- function(res,
                                spec_base,
                                spec_alt,
                                formula_base,
                                formula_alt,
                                comparison_name) {
    a <- get_spec_estimates(res, spec_base, formula_base)
    names(a)[names(a) == "estimate"] <- "base"
    b <- get_spec_estimates(res, spec_alt, formula_alt)
    names(b)[names(b) == "estimate"] <- "alt"

    z <- merge(a, b, by = c("scenario", "rep", "i", "j", "type", "true"))
    z <- z[complete.cases(z), , drop = FALSE]
    z <- z[
      is.finite(z$base) & is.finite(z$alt) & is.finite(z$true),
      , drop = FALSE
    ]
    if (nrow(z) == 0) return(NULL)

    z$comparison <- comparison_name
    z$delta <- z$alt - z$base
    z$abs_delta <- abs(z$delta)
    z$sq_delta <- z$delta^2

    key <- interaction(z$scenario, z$type, drop = TRUE)
    out <- lapply(split(z, key), function(d) {
      data.frame(
        scenario = d$scenario[1],
        comparison = comparison_name,
        type = d$type[1],
        n_obs = nrow(d),
        mean_delta = mean(d$delta),
        mean_abs_delta = mean(d$abs_delta),
        rmse_delta = sqrt(mean(d$sq_delta)),
        share_positive_delta = share_above_tolerance(d$delta)
      )
    })

    safe_rbind(out)
  }

  summarise_did <- function(res,
                            formula_la = "Green_Alston",
                            exact_base = "AIDS_P",
                            exact_uv = "AIDS_UV_norm",
                            la_base = "LA_PP_current",
                            la_uv = "LA_UVraw_UVraw_current",
                            label = "DID_raw_UV") {

    A0 <- get_spec_estimates(res, exact_base, "AIDS")
    names(A0)[names(A0) == "estimate"] <- "AIDS_base"

    A1 <- get_spec_estimates(res, exact_uv, "AIDS")
    names(A1)[names(A1) == "estimate"] <- "AIDS_uv"

    L0 <- get_spec_estimates(res, la_base, formula_la)
    names(L0)[names(L0) == "estimate"] <- "LA_base"

    L1 <- get_spec_estimates(res, la_uv, formula_la)
    names(L1)[names(L1) == "estimate"] <- "LA_uv"

    z <- merge(A0, A1, by = c("scenario", "rep", "i", "j", "type", "true"))
    z <- merge(z, L0, by = c("scenario", "rep", "i", "j", "type", "true"))
    z <- merge(z, L1, by = c("scenario", "rep", "i", "j", "type", "true"))
    z <- z[complete.cases(z), , drop = FALSE]
    finite_columns <- c("true", "AIDS_base", "AIDS_uv", "LA_base", "LA_uv")
    z <- z[
      apply(z[, finite_columns, drop = FALSE], 1, function(row) {
        all(is.finite(as.numeric(row)))
      }),
      , drop = FALSE
    ]
    if (nrow(z) == 0) return(NULL)

    z$comparison <- label
    z$delta_AIDS <- z$AIDS_uv - z$AIDS_base
    z$delta_LA <- z$LA_uv - z$LA_base
    z$extra_abs <- abs(z$delta_LA) - abs(z$delta_AIDS)
    z$extra_sq <- z$delta_LA^2 - z$delta_AIDS^2

    key <- interaction(z$scenario, z$type, drop = TRUE)
    out <- lapply(split(z, key), function(d) {
      data.frame(
        scenario = d$scenario[1],
        comparison = label,
        formula = formula_la,
        type = d$type[1],
        n_obs = nrow(d),
        mean_abs_delta_AIDS = mean(abs(d$delta_AIDS)),
        mean_abs_delta_LA = mean(abs(d$delta_LA)),
        mean_extra_abs = mean(d$extra_abs),
        mean_extra_sq = mean(d$extra_sq),
        share_LA_larger = share_above_tolerance(d$extra_abs)
      )
    })

    safe_rbind(out)
  }

  summarise_beta <- function(beta_df) {
    beta_df <- beta_df[
      is.finite(beta_df$beta_hat) & is.finite(beta_df$beta_true) &
        is.finite(beta_df$beta_error),
      , drop = FALSE
    ]
    if (nrow(beta_df) == 0) return(NULL)
    beta_df$abs_error <- abs(beta_df$beta_error)
    beta_df$sq_error <- beta_df$beta_error^2

    key <- interaction(
      beta_df$scenario,
      beta_df$spec,
      beta_df$model,
      beta_df$price_regressor,
      beta_df$index_price,
      beta_df$index_type,
      beta_df$i,
      drop = TRUE
    )

    out <- lapply(split(beta_df, key), function(z) {
      data.frame(
        scenario = z$scenario[1],
        spec = z$spec[1],
        model = z$model[1],
        price_regressor = z$price_regressor[1],
        index_price = z$index_price[1],
        index_type = z$index_type[1],
        i = z$i[1],
        beta_true = z$beta_true[1],
        mean_beta_hat = mean(z$beta_hat, na.rm = TRUE),
        beta_bias = mean(z$beta_error, na.rm = TRUE),
        mean_abs_beta_error = mean(z$abs_error, na.rm = TRUE),
        beta_rmse = sqrt(mean(z$sq_error, na.rm = TRUE))
      )
    })

    safe_rbind(out)
  }

  summarise_theoretical <- function(th_df) {
    key <- interaction(th_df$scenario, th_df$i, drop = TRUE)
    out <- lapply(split(th_df, key), function(z) {
      data.frame(
        scenario = z$scenario[1],
        i = z$i[1],
        n_rep = nrow(z),
        mean_abs_C_centered = mean(z$mean_abs_C_centered, na.rm = TRUE),
        sd_C_centered = mean(z$sd_C_centered, na.rm = TRUE),
        mean_abs_K_index = mean(z$mean_abs_K_index, na.rm = TRUE),
        mean_abs_tau_inner = mean(z$mean_abs_tau_inner, na.rm = TRUE),
        mean_abs_D_comm = mean(z$mean_abs_D_comm, na.rm = TRUE),
        mean_abs_D_tau = mean(z$mean_abs_D_tau, na.rm = TRUE),
        mean_abs_D_total = mean(z$mean_abs_D_total, na.rm = TRUE),
        rmse_D_comm = mean(z$rmse_D_comm, na.rm = TRUE),
        rmse_D_tau = mean(z$rmse_D_tau, na.rm = TRUE),
        rmse_D_total = mean(z$rmse_D_total, na.rm = TRUE),
        mcse_mean_abs_D_comm = sd(z$mean_abs_D_comm, na.rm = TRUE) / sqrt(sum(is.finite(z$mean_abs_D_comm))),
        mcse_mean_abs_D_tau = sd(z$mean_abs_D_tau, na.rm = TRUE) / sqrt(sum(is.finite(z$mean_abs_D_tau)))
      )
    })
    safe_rbind(out)
  }

  summarise_convergence <- function(conv_df) {
    key <- interaction(
      conv_df$scenario,
      conv_df$spec,
      conv_df$model,
      conv_df$price_regressor,
      conv_df$index_price,
      conv_df$index_type,
      drop = TRUE
    )

    out <- lapply(split(conv_df, key), function(z) {
      data.frame(
        scenario = z$scenario[1],
        spec = z$spec[1],
        model = z$model[1],
        price_regressor = z$price_regressor[1],
        index_price = z$index_price[1],
        index_type = z$index_type[1],
        n_rep = nrow(z),
        convergence_rate = mean(z$ok, na.rm = TRUE),
        mean_objective = if (any(is.finite(z$objective))) {
          mean(z$objective[is.finite(z$objective)])
        } else {
          NA_real_
        },
        median_objective = if (any(is.finite(z$objective))) {
          median(z$objective[is.finite(z$objective)])
        } else {
          NA_real_
        }
      )
    })

    safe_rbind(out)
  }

  # Scenario design
  make_scenarios <- function() {

    base_cfg <- list(
      seed = 12345,
      R = config$mc_reps,
      N = config$mc_households,

      beta_true = c(0.04, -0.02, 0.03, -0.05),
      gamma_scale = 1,
      homothetic = FALSE,

      sigma_logp = 0.20,
      price_corr = 0.30,

      M_mean = 0.00,
      M_sd = 0.40,

      sigma_u = 0.010,

      min_share = 0.02,
      max_tries = 200,

      quality_income = c(0.08, 0.04, -0.03, -0.09),
      quality_price  = c(0.20, 0.10, 0.15, 0.05),
      taste_loading  = 0.30,

      sigma_exp_error = 0.020,
      sigma_qty_error = 0.040,
      sigma_extra_uv  = 0.020,
      share_reporting_error = TRUE,

      n_clusters = 60,
      sigma_local_cluster = 0.000,
      sigma_local_household = 0.000,

      unit_scale = rep(0, 4),

      run_invariant_indices = TRUE,

      verbose = TRUE,
      print_every = 50
    )

    moderate_scale <- c(log(10), log(5), log(2), 0)
    large_scale <- c(log(1000), log(100), log(10), 0)
    mixed_scale <- c(log(1000), -log(10), log(50), 0)

    scenarios <- list(
      normalized_baseline = base_cfg,

      unit_scale_moderate = modifyList(
        base_cfg,
        list(
          seed = 22345,
          unit_scale = moderate_scale,
          sigma_exp_error = 0.010,
          sigma_qty_error = 0.020,
          sigma_extra_uv = 0.010,
          sigma_local_cluster = 0.000,
          sigma_local_household = 0.000
        )
      ),

      unit_scale_large = modifyList(
        base_cfg,
        list(
          seed = 32345,
          unit_scale = large_scale,
          sigma_exp_error = 0.010,
          sigma_qty_error = 0.020,
          sigma_extra_uv = 0.010,
          sigma_local_cluster = 0.000,
          sigma_local_household = 0.000
        )
      ),

      unit_scale_mixed = modifyList(
        base_cfg,
        list(
          seed = 42345,
          unit_scale = mixed_scale,
          sigma_exp_error = 0.010,
          sigma_qty_error = 0.020,
          sigma_extra_uv = 0.010,
          sigma_local_cluster = 0.000,
          sigma_local_household = 0.000
        )
      ),

      local_unit_error = modifyList(
        base_cfg,
        list(
          seed = 52345,
          unit_scale = rep(0, 4),
          sigma_local_cluster = 0.120,
          sigma_local_household = 0.040,
          sigma_exp_error = 0.010,
          sigma_qty_error = 0.020,
          sigma_extra_uv = 0.010
        )
      ),

      quality_reporting_only = modifyList(
        base_cfg,
        list(
          seed = 62345,
          unit_scale = rep(0, 4),
          quality_income = c(0.16, 0.08, -0.06, -0.18),
          quality_price = c(0.35, 0.20, 0.25, 0.10),
          taste_loading = 0.50,
          sigma_exp_error = 0.050,
          sigma_qty_error = 0.080,
          sigma_extra_uv = 0.040,
          sigma_local_cluster = 0.000,
          sigma_local_household = 0.000
        )
      ),

      combined_unit_quality = modifyList(
        base_cfg,
        list(
          seed = 72345,
          unit_scale = moderate_scale,
          quality_income = c(0.16, 0.08, -0.06, -0.18),
          quality_price = c(0.35, 0.20, 0.25, 0.10),
          taste_loading = 0.50,
          sigma_exp_error = 0.050,
          sigma_qty_error = 0.080,
          sigma_extra_uv = 0.040,
          sigma_local_cluster = 0.120,
          sigma_local_household = 0.040
        )
      ),

      high_beta_unit_scale = modifyList(
        base_cfg,
        list(
          seed = 82345,
          beta_true = c(0.10, -0.06, 0.08, -0.12),
          unit_scale = moderate_scale,
          M_sd = 0.30,
          sigma_u = 0.008,
          sigma_exp_error = 0.010,
          sigma_qty_error = 0.020,
          sigma_extra_uv = 0.010
        )
      ),

      fixed_scale_only_mixed = modifyList(
        base_cfg,
        list(
          seed = 92345,
          unit_scale = mixed_scale,
          quality_income = rep(0, 4),
          quality_price = rep(0, 4),
          taste_loading = 0,
          sigma_exp_error = 0,
          sigma_qty_error = 0,
          sigma_extra_uv = 0,
          share_reporting_error = FALSE,
          sigma_local_cluster = 0,
          sigma_local_household = 0
        )
      ),

      homothetic_fixed_scale_mixed = modifyList(
        base_cfg,
        list(
          seed = 102345,
          beta_true = rep(0, 4),
          homothetic = TRUE,
          # Disable demand noise as well as beta for the homothetic null.
          sigma_u = 0,
          unit_scale = mixed_scale,
          quality_income = rep(0, 4),
          quality_price = rep(0, 4),
          taste_loading = 0,
          sigma_exp_error = 0,
          sigma_qty_error = 0,
          sigma_extra_uv = 0,
          share_reporting_error = FALSE,
          sigma_local_cluster = 0,
          sigma_local_household = 0
        )
      )
    )

    scenarios
  }

  # Monte Carlo execution
  canonical_scenarios <- make_scenarios()
  requested_scenario_filter <- config$mc_scenario_filter

  if (length(requested_scenario_filter) > 0) {
    unknown_scenarios <- setdiff(requested_scenario_filter, names(canonical_scenarios))
    if (length(unknown_scenarios) > 0) {
      stop(
        "Unknown mc_scenario_filter value(s): ",
        paste(unknown_scenarios, collapse = ", "),
        "\nAvailable scenarios: ",
        paste(names(canonical_scenarios), collapse = ", "),
        call. = FALSE
      )
    }
    scenarios <- canonical_scenarios[requested_scenario_filter]
    message(
      "Monte Carlo scenario filter active: ",
      paste(names(scenarios), collapse = ", ")
    )
  } else {
    scenarios <- canonical_scenarios
  }

  all_results <- vector("list", length(scenarios))
  names(all_results) <- names(scenarios)
  for (nm in names(scenarios)) {
    message("Monte Carlo: ", nm)
    all_results[[nm]] <- run_monte_carlo_scenario(nm, scenarios[[nm]])
  }
  elasticities_all <- safe_rbind(lapply(all_results, `[[`, "elasticities"))
  beta_all <- safe_rbind(lapply(all_results, `[[`, "beta"))
  convergence_all <- safe_rbind(lapply(all_results, `[[`, "convergence"))
  theoretical_all <- safe_rbind(lapply(all_results, `[[`, "theoretical"))
  oracle_all <- safe_rbind(lapply(all_results, `[[`, "oracle"))
  replication_execution_all <- safe_rbind(lapply(all_results, `[[`, "replication_execution"))
  true_parameters_by_scenario <- lapply(all_results, `[[`, "true_parameters")
  scenarios <- lapply(all_results, `[[`, "configuration")
  mc_raw_created_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  raw_results <- list(
    format_version = "1.0.0",
    created_at = mc_raw_created_at,
    rng_kind = RNGkind(),
    configuration_provenance = mc_configuration_provenance,
    scenario_configurations = scenarios,
    true_parameters_by_scenario = true_parameters_by_scenario,
    replication_execution = replication_execution_all,
    elasticities = elasticities_all,
    beta = beta_all,
    convergence = convergence_all,
    theoretical_by_replication = theoretical_all,
    oracle = oracle_all
  )
  saveRDS(raw_results, mc_raw_results_path, compress = TRUE)
  write.csv(replication_execution_all,
            file.path(output_dir, "mc_replication_execution_diagnostics.csv"), row.names = FALSE)
  if (is.null(elasticities_all) || nrow(elasticities_all) == 0 ||
      is.null(beta_all) || nrow(beta_all) == 0 ||
      is.null(convergence_all) || nrow(convergence_all) == 0) {
    stop("Monte Carlo produced no usable model records; see mc_replication_execution_diagnostics.csv.",
         call. = FALSE)
  }

  # Tables and diagnostics
  summary_elasticities <- summarise_elasticities(elasticities_all)
  summary_beta <- summarise_beta(beta_all)
  summary_convergence <- summarise_convergence(convergence_all)
  summary_theoretical <- summarise_theoretical(theoretical_all)

  did_general_GA <- summarise_did(
    elasticities_all,
    formula_la = "Green_Alston",
    exact_base = "AIDS_P",
    exact_uv = "AIDS_UV_norm",
    la_base = "LA_PP_current",
    la_uv = "LA_UVnorm_UVnorm_current",
    label = "DID_general_UV_normalized"
  )

  did_general_LAprime <- summarise_did(
    elasticities_all,
    formula_la = "LA_prime",
    exact_base = "AIDS_P",
    exact_uv = "AIDS_UV_norm",
    la_base = "LA_PP_current",
    la_uv = "LA_UVnorm_UVnorm_current",
    label = "DID_general_UV_normalized"
  )

  did_raw_GA <- summarise_did(
    elasticities_all,
    formula_la = "Green_Alston",
    exact_base = "AIDS_P",
    exact_uv = "AIDS_UV_norm",
    la_base = "LA_PP_current",
    la_uv = "LA_UVraw_UVraw_current",
    label = "DID_raw_UV_commensurability"
  )

  did_raw_LAprime <- summarise_did(
    elasticities_all,
    formula_la = "LA_prime",
    exact_base = "AIDS_P",
    exact_uv = "AIDS_UV_norm",
    la_base = "LA_PP_current",
    la_uv = "LA_UVraw_UVraw_current",
    label = "DID_raw_UV_commensurability"
  )

  stone_channel_raw_GA <- compare_two_specs(
    elasticities_all,
    spec_base = "LA_PP_current",
    spec_alt = "LA_P_UVraw_current",
    formula_base = "Green_Alston",
    formula_alt = "Green_Alston",
    comparison_name = "Stone_channel_only_raw_UV: LA(P,UV_raw)-LA(P,P)"
  )

  stone_channel_norm_GA <- compare_two_specs(
    elasticities_all,
    spec_base = "LA_PP_current",
    spec_alt = "LA_P_UVnorm_current",
    formula_base = "Green_Alston",
    formula_alt = "Green_Alston",
    comparison_name = "Stone_channel_only_norm_UV: LA(P,UV_norm)-LA(P,P)"
  )

  unit_scale_only_stone_GA <- compare_two_specs(
    elasticities_all,
    spec_base = "LA_P_UVnorm_current",
    spec_alt = "LA_P_UVraw_current",
    formula_base = "Green_Alston",
    formula_alt = "Green_Alston",
    comparison_name = "Unit_scale_only_in_Stone: LA(P,UV_raw)-LA(P,UV_norm)"
  )

  total_raw_LA_GA <- compare_two_specs(
    elasticities_all,
    spec_base = "LA_PP_current",
    spec_alt = "LA_UVraw_UVraw_current",
    formula_base = "Green_Alston",
    formula_alt = "Green_Alston",
    comparison_name = "Total_raw_UV_in_LA: LA(UV_raw,UV_raw)-LA(P,P)"
  )

  direct_raw_regressor_GA <- compare_two_specs(
    elasticities_all,
    spec_base = "LA_PP_current",
    spec_alt = "LA_UVraw_P_current",
    formula_base = "Green_Alston",
    formula_alt = "Green_Alston",
    comparison_name = "Direct_regressor_raw_UV_only: LA(UV_raw,P)-LA(P,P)"
  )

  stone_decomposition_GA <- rbind(
    stone_channel_raw_GA,
    stone_channel_norm_GA,
    unit_scale_only_stone_GA,
    total_raw_LA_GA,
    direct_raw_regressor_GA
  )

  invariant_comparisons_GA <- NULL
  if (any(elasticities_all$spec == "LA_UVraw_UVraw_corrected")) {
    invariant_comparisons_GA <- rbind(
      compare_two_specs(
        elasticities_all,
        spec_base = "LA_UVnorm_UVnorm_corrected",
        spec_alt = "LA_UVraw_UVraw_corrected",
        formula_base = "Green_Alston",
        formula_alt = "Green_Alston",
        comparison_name = "Corrected_index_unit_scale_sensitivity: raw-norm"
      ),
      compare_two_specs(
        elasticities_all,
        spec_base = "LA_UVnorm_UVnorm_tornqvist",
        spec_alt = "LA_UVraw_UVraw_tornqvist",
        formula_base = "Green_Alston",
        formula_alt = "Green_Alston",
        comparison_name = "Tornqvist_index_unit_scale_sensitivity: raw-norm"
      )
    )
  }

  paper_table_commensurability_GA <- rbind(
    did_general_GA,
    did_raw_GA
  )
  if (isTRUE(config$mc_write_raw_csv)) {
    write.csv(elasticities_all,
              file.path(output_dir, "mc_elasticities_raw.csv"),
              row.names = FALSE)
    write.csv(beta_all,
              file.path(output_dir, "mc_beta_raw.csv"),
              row.names = FALSE)
  }
  write.csv(convergence_all,
            file.path(output_dir, "mc_convergence_raw.csv"),
            row.names = FALSE)
  write.csv(theoretical_all,
            file.path(output_dir, "mc_theoretical_channels_by_replication.csv"),
            row.names = FALSE)
  write.csv(oracle_all,
            file.path(output_dir, "mc_oracle_exact_AIDS_commensurability_check.csv"),
            row.names = FALSE)

  write.csv(summary_elasticities,
            file.path(output_dir, "mc_summary_elasticities_by_scenario.csv"),
            row.names = FALSE)
  write.csv(summary_beta,
            file.path(output_dir, "mc_beta_summary_by_scenario.csv"),
            row.names = FALSE)
  write.csv(summary_convergence,
            file.path(output_dir, "mc_convergence_summary_by_scenario.csv"),
            row.names = FALSE)
  write.csv(summary_theoretical,
            file.path(output_dir, "mc_theoretical_channels_summary_by_scenario.csv"),
            row.names = FALSE)

  write.csv(did_general_GA,
            file.path(output_dir, "mc_did_general_uv_green_alston_by_scenario.csv"),
            row.names = FALSE)
  write.csv(did_general_LAprime,
            file.path(output_dir, "mc_did_general_uv_la_prime_by_scenario.csv"),
            row.names = FALSE)
  write.csv(did_raw_GA,
            file.path(output_dir, "mc_did_raw_uv_commensurability_green_alston_by_scenario.csv"),
            row.names = FALSE)
  write.csv(did_raw_LAprime,
            file.path(output_dir, "mc_did_raw_uv_commensurability_la_prime_by_scenario.csv"),
            row.names = FALSE)

  write.csv(stone_decomposition_GA,
            file.path(output_dir, "mc_la_stone_channel_decomposition_green_alston_by_scenario.csv"),
            row.names = FALSE)

  if (!is.null(invariant_comparisons_GA)) {
    write.csv(invariant_comparisons_GA,
              file.path(output_dir, "mc_invariant_index_unit_scale_sensitivity_green_alston.csv"),
              row.names = FALSE)
  }

  write.csv(paper_table_commensurability_GA,
            file.path(output_dir, "paper_table_commensurability_green_alston.csv"),
            row.names = FALSE)

  replication_level_spec_difference <- function(res,
                                                spec_base,
                                                spec_alt,
                                                formula_name,
                                                label) {
    base <- subset(
      res,
      spec == spec_base & formula == formula_name,
      select = c("scenario", "rep", "i", "j", "type", "estimate")
    )
    names(base)[names(base) == "estimate"] <- "base"
    alt <- subset(
      res,
      spec == spec_alt & formula == formula_name,
      select = c("scenario", "rep", "i", "j", "type", "estimate")
    )
    names(alt)[names(alt) == "estimate"] <- "alt"
    z <- merge(base, alt, by = c("scenario", "rep", "i", "j", "type"))
    z <- z[complete.cases(z), ]
    z$abs_difference <- abs(z$alt - z$base)

    key <- interaction(z$scenario, z$rep, z$type, drop = TRUE)
    per_rep <- lapply(split(z, key), function(d) {
      data.frame(
        scenario = d$scenario[1],
        rep = d$rep[1],
        comparison = label,
        formula = formula_name,
        type = d$type[1],
        mean_abs_difference = mean(d$abs_difference)
      )
    })
    safe_rbind(per_rep)
  }

  summarise_replication_level <- function(rep_df) {
    key <- interaction(
      rep_df$scenario, rep_df$comparison, rep_df$formula, rep_df$type,
      drop = TRUE
    )
    out <- lapply(split(rep_df, key), function(d) {
      n_rep <- nrow(d)
      data.frame(
        scenario = d$scenario[1],
        comparison = d$comparison[1],
        formula = d$formula[1],
        type = d$type[1],
        n_rep = n_rep,
        mean = mean(d$mean_abs_difference),
        sd_across_replications = sd(d$mean_abs_difference),
        monte_carlo_se = sd(d$mean_abs_difference) / sqrt(n_rep),
        p025 = as.numeric(quantile(d$mean_abs_difference, 0.025, names = FALSE)),
        p975 = as.numeric(quantile(d$mean_abs_difference, 0.975, names = FALSE))
      )
    })
    safe_rbind(out)
  }

  unit_scale_stone_rep <- replication_level_spec_difference(
    elasticities_all,
    spec_base = "LA_P_UVnorm_current",
    spec_alt = "LA_P_UVraw_current",
    formula_name = "Green_Alston",
    label = "Unit_scale_only_in_Stone"
  )
  current_raw_total_rep <- replication_level_spec_difference(
    elasticities_all,
    spec_base = "LA_PP_current",
    spec_alt = "LA_UVraw_UVraw_current",
    formula_name = "Green_Alston",
    label = "Total_raw_UV_in_LA"
  )
  replication_level_key_comparisons <- rbind(unit_scale_stone_rep, current_raw_total_rep)
  replication_level_key_summary <- summarise_replication_level(replication_level_key_comparisons)

  beta_replication_summary <- NULL
  if (!is.null(beta_all) && nrow(beta_all) > 0) {
    beta_ok <- beta_all[is.finite(beta_all$beta_hat), ]
    key <- interaction(beta_ok$scenario, beta_ok$spec, beta_ok$i, drop = TRUE)
    beta_replication_summary <- safe_rbind(lapply(split(beta_ok, key), function(d) {
      n_rep <- nrow(d)
      data.frame(
        scenario = d$scenario[1],
        spec = d$spec[1],
        i = d$i[1],
        beta_true = d$beta_true[1],
        n_rep = n_rep,
        mean_beta_hat = mean(d$beta_hat),
        sd_beta_hat = sd(d$beta_hat),
        monte_carlo_se_beta_hat = sd(d$beta_hat) / sqrt(n_rep),
        p025_beta_hat = as.numeric(quantile(d$beta_hat, 0.025, names = FALSE)),
        p975_beta_hat = as.numeric(quantile(d$beta_hat, 0.975, names = FALSE))
      )
    }))
  }

  # Design and execution metadata
  scenario_purpose <- c(
    normalized_baseline = "General unit-value contamination without fixed commodity unit scales.",
    unit_scale_moderate = "Moderate fixed unit scales with small residual unit-value noise.",
    unit_scale_large = "Large fixed unit scales with small residual unit-value noise.",
    unit_scale_mixed = "Mixed-sign fixed unit scales; stresses relative-unit sensitivity.",
    local_unit_error = "Household and cluster conversion error without fixed commodity scales.",
    quality_reporting_only = "Stronger quality and reporting components without fixed commodity scales.",
    combined_unit_quality = "Fixed scales combined with quality, reporting, and local conversion components.",
    high_beta_unit_scale = "Calibrated larger-beta sensitivity design, also using smaller expenditure and demand-disturbance variances.",
    fixed_scale_only_mixed = "Clean mechanism scenario: fixed scales only, with all tau components set to zero.",
    homothetic_fixed_scale_mixed = "Deterministic homothetic null/sanity check: fixed scales are present, beta equals zero, and the demand disturbance is disabled; the Stone real-expenditure channel should vanish up to numerical tolerance."
  )

  mc_format_value <- function(x) {
    if (is.null(x)) return("NULL")
    if (is.numeric(x)) return(paste(sprintf("%.17g", as.numeric(x)), collapse = ";"))
    if (is.logical(x)) return(paste(ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE")), collapse = ";"))
    paste(as.character(x), collapse = ";")
  }

  mc_parameter_groups <- c(
    seed = "execution",
    R = "execution",
    N = "execution",
    verbose = "execution",
    print_every = "execution",
    run_invariant_indices = "estimation_design",
    beta_true = "true_AIDS_parameters",
    gamma_scale = "true_AIDS_parameters",
    homothetic = "true_AIDS_parameters",
    sigma_logp = "standardized_price_process",
    price_corr = "standardized_price_process",
    M_mean = "real_expenditure_process",
    M_sd = "real_expenditure_process",
    sigma_u = "demand_disturbance",
    min_share = "share_admissibility",
    max_tries = "share_admissibility",
    quality_income = "unit_value_quality_component",
    quality_price = "unit_value_quality_component",
    taste_loading = "unit_value_quality_component",
    sigma_exp_error = "reporting_error_component",
    sigma_qty_error = "reporting_error_component",
    sigma_extra_uv = "other_unit_value_noise",
    share_reporting_error = "reporting_error_component",
    n_clusters = "local_conversion_component",
    sigma_local_cluster = "local_conversion_component",
    sigma_local_household = "local_conversion_component",
    unit_scale = "fixed_measurement_unit_scale"
  )

  mc_scenario_parameter_manifest <- safe_rbind(lapply(names(scenarios), function(nm) {
    cfg <- scenarios[[nm]]
    safe_rbind(lapply(names(cfg), function(parameter_name) {
      value <- cfg[[parameter_name]]
      data.frame(
        scenario = nm,
        parameter_group = if (parameter_name %in% names(mc_parameter_groups)) {
          unname(mc_parameter_groups[parameter_name])
        } else {
          "other"
        },
        parameter = parameter_name,
        value_class = class(value)[1],
        vector_length = length(value),
        value = mc_format_value(value),
        configuration_provenance = mc_configuration_provenance,
        stringsAsFactors = FALSE
      )
    }))
  }))

  mc_model_catalog <- data.frame(
    spec = c(
      "AIDS_P", "AIDS_UV_norm", "LA_PP_current",
      "LA_UVnorm_UVnorm_current", "LA_UVraw_UVraw_current",
      "LA_P_UVnorm_current", "LA_P_UVraw_current",
      "LA_UVnorm_P_current", "LA_UVraw_P_current",
      "LA_UVraw_UVraw_corrected", "LA_UVraw_UVraw_tornqvist",
      "LA_UVnorm_UVnorm_corrected", "LA_UVnorm_UVnorm_tornqvist"
    ),
    model = c(
      "AIDS", "AIDS", rep("LA_AIDS", 11)
    ),
    price_regressor = c(
      "P", "UV_norm", "P", "UV_norm", "UV_raw", "P", "P",
      "UV_norm", "UV_raw", "UV_raw", "UV_raw", "UV_norm", "UV_norm"
    ),
    index_price = c(
      "translog", "translog", "P", "UV_norm", "UV_raw", "UV_norm",
      "UV_raw", "P", "P", "UV_raw", "UV_raw", "UV_norm", "UV_norm"
    ),
    index_type = c(
      "translog", "translog", rep("stone_current", 7),
      "corrected_stone", "tornqvist", "corrected_stone", "tornqvist"
    ),
    exact_model = c(TRUE, TRUE, rep(FALSE, 11)),
    optional_invariant_index_specification = c(rep(FALSE, 9), rep(TRUE, 4)),
    elasticity_formulas_stored = c(
      "AIDS", "AIDS", rep("Green_Alston;LA_prime", 11)
    ),
    stringsAsFactors = FALSE
  )

  mc_model_specification_manifest <- safe_rbind(lapply(names(scenarios), function(nm) {
    cfg <- scenarios[[nm]]
    keep <- !mc_model_catalog$optional_invariant_index_specification |
      isTRUE(cfg$run_invariant_indices)
    z <- mc_model_catalog[keep, , drop = FALSE]
    z$scenario <- nm
    z$objective <- ifelse(
      z$exact_model,
      "Restricted nonlinear least squares on the first n-1 share equations using the exact translog price aggregator.",
      "Restricted nonlinear least squares on the first n-1 share equations using the stated linear price index."
    )
    z$parameterization <- paste(
      "alpha0 fixed at zero; adding-up, homogeneity, and symmetry imposed",
      "through alpha, beta, and the symmetric upper-left Gamma block"
    )
    z$starting_values <- "Equation-by-equation LA starting values, symmetrized Gamma block, followed by deterministic retry starts."
    z$convergence_rule <- "optim convergence code 0; non-converged estimates retained only in diagnostics and excluded from numerical summaries."
    z$configuration_provenance <- mc_configuration_provenance
    z[, c(
      "scenario", "spec", "model", "exact_model", "price_regressor",
      "index_price", "index_type", "elasticity_formulas_stored",
      "optional_invariant_index_specification", "objective",
      "parameterization", "starting_values", "convergence_rule",
      "configuration_provenance"
    )]
  }))

  mc_optimizer_manifest <- data.frame(
    setting = c(
      "estimator", "objective_equations", "alpha0_treatment",
      "restriction_parameterization", "starting_values",
      "attempt_1", "attempt_2", "attempt_3", "attempt_4",
      "deterministic_retry_shift_scale", "maximum_iterations",
      "relative_tolerance", "convergence_acceptance_rule",
      "nonconverged_estimate_treatment", "evaluation_shares",
      "evaluation_prices", "random_number_stream"
    ),
    value = c(
      "Restricted nonlinear least squares implemented with optim",
      "First n-1 budget-share equations; the omitted equation follows by adding-up",
      "Fixed at zero after price normalization in every estimated exact/linear model",
      "Free first n-1 alpha and beta elements plus the symmetric upper-left Gamma block; adding-up, homogeneity, and symmetry completed algebraically",
      "Equation-by-equation linear approximation, with the Gamma block symmetrized",
      "BFGS from the starting vector",
      "BFGS from starting vector plus deterministic shift",
      "BFGS from starting vector minus deterministic shift",
      "Nelder-Mead from the original starting vector",
      format(mc_retry_shift_scale, scientific = TRUE),
      as.character(mc_optimizer_maxit),
      format(mc_optimizer_reltol, scientific = TRUE),
      "optim convergence code exactly zero",
      "Stored in convergence diagnostics but excluded from elasticity and coefficient summaries",
      "Column means of observed budget shares in each replication",
      "Column means of the specification-specific price regressors/index prices in each replication",
      "One scenario-level seed followed by a sequential R RNG stream within the scenario"
    ),
    code_function = c(
      "fit_aids_model", "fit_aids_model", "unpack_theta",
      "pack_theta;unpack_theta;complete_gamma_from_block", "make_start_from_la",
      rep("fit_aids_model", 9),
      "run_one_replication", "run_one_replication", "run_monte_carlo_scenario"
    ),
    configuration_provenance = mc_configuration_provenance,
    stringsAsFactors = FALSE
  )

  mc_true_parameter_manifest <- safe_rbind(lapply(names(scenarios), function(nm) {
    par <- true_parameters_by_scenario[[nm]]
    if (is.null(par)) return(NULL)
    rows <- list(
      data.frame(
        scenario = nm, parameter = "n", i = NA_integer_, j = NA_integer_,
        value = as.numeric(par$n), parameter_provenance = mc_true_parameter_provenance
      ),
      data.frame(
        scenario = nm, parameter = "alpha0", i = NA_integer_, j = NA_integer_,
        value = as.numeric(par$alpha0), parameter_provenance = mc_true_parameter_provenance
      ),
      data.frame(
        scenario = nm, parameter = "alpha", i = seq_along(par$alpha),
        j = NA_integer_, value = as.numeric(par$alpha),
        parameter_provenance = mc_true_parameter_provenance
      ),
      data.frame(
        scenario = nm, parameter = "beta", i = seq_along(par$beta),
        j = NA_integer_, value = as.numeric(par$beta),
        parameter_provenance = mc_true_parameter_provenance
      )
    )
    gamma_grid <- expand.grid(
      i = seq_len(nrow(par$Gamma)),
      j = seq_len(ncol(par$Gamma))
    )
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = nm,
      parameter = "Gamma",
      i = gamma_grid$i,
      j = gamma_grid$j,
      value = as.vector(par$Gamma),
      parameter_provenance = mc_true_parameter_provenance
    )
    safe_rbind(rows)
  }))

  mc_dgp_equation_manifest <- data.frame(
    step = c(
      "standardized_log_prices", "true_real_expenditure", "demand_disturbance",
      "true_budget_shares", "true_log_expenditure", "income_quality",
      "price_quality", "taste_quality", "reporting_component",
      "local_conversion_component", "other_unit_value_noise",
      "normalized_unit_values", "raw_unit_values", "observed_budget_shares",
      "observed_log_expenditure", "Stone_fixed_scale_component"
    ),
    object = c(
      "L_true", "M", "U", "W_true", "m_true", "Tau_income", "Tau_price",
      "Tau_taste", "Tau_reporting", "Kappa", "Tau_extra", "L_uv_norm",
      "L_uv_raw", "W_obs", "m_obs", "C_h_obs"
    ),
    definition = c(
      "Rows are mean-zero multivariate normal with variance sigma_logp^2 and equicorrelation price_corr.",
      "M_h is normal with mean M_mean and standard deviation M_sd.",
      "Independent normal draws with standard deviation sigma_u are demeaned across goods so each row sums to zero.",
      "alpha + Gamma L_true + beta M + U; data are redrawn until every share lies between min_share and 1-min_share.",
      "M + a(L_true), where a(.) is the exact translog AIDS price aggregator with alpha0 fixed at zero.",
      "Outer product of M and quality_income.",
      "Minus L_true multiplied good by good by quality_price.",
      "taste_loading multiplied by U.",
      "A_exp - B_qty, with independent normal expenditure and quantity errors.",
      "Cluster-specific normal component plus household-specific normal component; households are assigned to n_clusters with replacement.",
      "Independent normal residual unit-value noise with standard deviation sigma_extra_uv.",
      "L_true + Tau; the fixed commodity scale vector is removed but household-specific contamination remains.",
      "L_true + unit_scale + Tau.",
      "If share_reporting_error is TRUE, W_true*exp(A_exp) normalized to sum to one; otherwise W_true.",
      "m_true plus the log of the expenditure-reporting normalization factor.",
      "W_obs multiplied by the fixed unit_scale vector."
    ),
    code_function = c(
      rep("generate_data", 16)
    ),
    stringsAsFactors = FALSE
  )

  actual_replications <- vapply(
    names(scenarios),
    function(nm) length(unique(convergence_all$rep[as.character(convergence_all$scenario) == nm])),
    integer(1)
  )

  scenario_manifest <- safe_rbind(lapply(names(scenarios), function(nm) {
    cfg <- scenarios[[nm]]
    exec <- replication_execution_all[
      as.character(replication_execution_all$scenario) == nm,
      , drop = FALSE
    ]
    data.frame(
      scenario = nm,
      purpose = unname(scenario_purpose[nm]),
      replications = actual_replications[nm],
      replications_requested = as.integer(cfg$R),
      replications_with_model_records = actual_replications[nm],
      replications_all_models_converged = sum(
        grepl("all_models_converged", exec$status), na.rm = TRUE
      ),
      complete_replication_failures_recorded = sum(
        exec$status == "failed_before_model_records", na.rm = TRUE
      ),
      households_per_replication = as.integer(cfg$N),
      seed = as.integer(cfg$seed),
      alpha_vector = mc_format_value(true_parameters_by_scenario[[nm]]$alpha),
      beta_vector = mc_format_value(cfg$beta_true),
      gamma_scale = as.numeric(cfg$gamma_scale),
      unit_scale_vector = mc_format_value(cfg$unit_scale),
      fixed_unit_scale_present = any(abs(cfg$unit_scale) > 1e-12),
      homothetic = isTRUE(cfg$homothetic),
      sigma_logp = as.numeric(cfg$sigma_logp),
      price_corr = as.numeric(cfg$price_corr),
      M_mean = as.numeric(cfg$M_mean),
      M_sd = as.numeric(cfg$M_sd),
      sigma_u = as.numeric(cfg$sigma_u),
      min_share = as.numeric(cfg$min_share),
      max_tries = as.integer(cfg$max_tries),
      quality_income_vector = mc_format_value(cfg$quality_income),
      quality_price_vector = mc_format_value(cfg$quality_price),
      taste_loading = as.numeric(cfg$taste_loading),
      sigma_exp_error = as.numeric(cfg$sigma_exp_error),
      sigma_qty_error = as.numeric(cfg$sigma_qty_error),
      sigma_extra_uv = as.numeric(cfg$sigma_extra_uv),
      share_reporting_error = isTRUE(cfg$share_reporting_error),
      n_clusters = as.integer(cfg$n_clusters),
      sigma_local_cluster = as.numeric(cfg$sigma_local_cluster),
      sigma_local_household = as.numeric(cfg$sigma_local_household),
      run_invariant_indices = isTRUE(cfg$run_invariant_indices),
      model_specification_count = sum(
        as.character(mc_model_specification_manifest$scenario) == nm
      ),
      numerical_zero_tolerance = mc_numerical_tol,
      configuration_provenance = mc_configuration_provenance,
      stringsAsFactors = FALSE
    )
  }))

  mc_replication_execution_summary <- safe_rbind(lapply(
    split(replication_execution_all, replication_execution_all$scenario),
    function(z) {
      status_table <- as.data.frame(table(z$status), stringsAsFactors = FALSE)
      names(status_table) <- c("status", "n_replications")
      status_table$scenario <- as.character(z$scenario[1])
      status_table[, c("scenario", "status", "n_replications")]
    }
  ))

  mc_reproducibility_manifest <- data.frame(
    setting = c("format_version", "created_at", "raw_result_path", "rng_kind",
                "numerical_zero_tolerance", "scenarios", "replications_per_scenario",
                "households_per_replication"),
    value = c("1.0.0", mc_raw_created_at,
              normalizePath(mc_raw_results_path, winslash = "/", mustWork = FALSE),
              mc_rng_kind_recorded, format(mc_numerical_tol, scientific = TRUE),
              paste(names(scenarios), collapse = ";"), as.character(config$mc_reps),
              as.character(config$mc_households)),
    stringsAsFactors = FALSE
  )

  # Replication uncertainty and Figure 1
  write.csv(
    replication_level_key_comparisons,
    file.path(output_dir, "mc_replication_level_key_comparisons.csv"),
    row.names = FALSE
  )
  write.csv(
    replication_level_key_summary,
    file.path(output_dir, "mc_replication_level_key_comparisons_with_mcse.csv"),
    row.names = FALSE
  )
  if (!is.null(beta_replication_summary)) {
    write.csv(
      beta_replication_summary,
      file.path(output_dir, "mc_beta_summary_with_mcse.csv"),
      row.names = FALSE
    )
  }
  write.csv(
    scenario_manifest,
    file.path(output_dir, "mc_scenario_manifest.csv"),
    row.names = FALSE
  )
  write.csv(
    mc_scenario_parameter_manifest,
    file.path(output_dir, "mc_complete_scenario_parameter_manifest.csv"),
    row.names = FALSE
  )
  write.csv(
    mc_true_parameter_manifest,
    file.path(output_dir, "mc_true_AIDS_parameter_manifest.csv"),
    row.names = FALSE
  )
  write.csv(
    mc_model_specification_manifest,
    file.path(output_dir, "mc_model_specification_manifest.csv"),
    row.names = FALSE
  )
  write.csv(
    mc_optimizer_manifest,
    file.path(output_dir, "mc_optimizer_and_numerical_method_manifest.csv"),
    row.names = FALSE
  )
  write.csv(
    mc_dgp_equation_manifest,
    file.path(output_dir, "mc_DGP_equation_manifest.csv"),
    row.names = FALSE
  )
  write.csv(
    mc_replication_execution_summary,
    file.path(output_dir, "mc_replication_execution_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    mc_reproducibility_manifest,
    file.path(output_dir, "mc_reproducibility_manifest.csv"),
    row.names = FALSE
  )

  falsification_audit <- subset(
    replication_level_key_summary,
    scenario == "homothetic_fixed_scale_mixed" &
      comparison == "Unit_scale_only_in_Stone"
  )
  if (nrow(falsification_audit) > 0) {
    falsification_tolerance <- 1e-6
    falsification_audit$theoretical_beta <- 0
    falsification_audit$demand_disturbance_sd <- unname(
      vapply(
        scenarios["homothetic_fixed_scale_mixed"],
        `[[`, numeric(1), "sigma_u"
      )
    )
    falsification_audit$tolerance <- falsification_tolerance
    falsification_audit$pass <- falsification_audit$mean <= falsification_tolerance
    write.csv(
      falsification_audit,
      file.path(output_dir, "mc_homothetic_falsification_audit.csv"),
      row.names = FALSE
    )
    if (any(!falsification_audit$pass)) {
      warning(
        "The deterministic homothetic null check exceeded the tolerance. ",
        "Inspect mc_homothetic_falsification_audit.csv before using this scenario."
      )
    }
  }

  make_tolerance_indicator_audit <- function(df, source_table, indicator_col) {
    z <- subset(df, scenario == "homothetic_fixed_scale_mixed")
    if (nrow(z) == 0) {
      return(data.frame(
        source_table = character(0),
        scenario = character(0),
        comparison = character(0),
        type = character(0),
        indicator = character(0),
        indicator_value = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      source_table = source_table,
      scenario = z$scenario,
      comparison = z$comparison,
      type = z$type,
      indicator = indicator_col,
      indicator_value = z[[indicator_col]],
      stringsAsFactors = FALSE
    )
  }

  tolerance_indicator_audit <- rbind(
    make_tolerance_indicator_audit(
      paper_table_commensurability_GA,
      "paper_table_commensurability_green_alston",
      "share_LA_larger"
    ),
    make_tolerance_indicator_audit(
      stone_decomposition_GA,
      "mc_la_stone_channel_decomposition_green_alston_by_scenario",
      "share_positive_delta"
    )
  )
  if (nrow(tolerance_indicator_audit) > 0) {
    tolerance_indicator_audit$numerical_tolerance <- mc_numerical_tol
    tolerance_indicator_audit$pass_zero_for_homothetic <-
      tolerance_indicator_audit$indicator_value == 0
    write.csv(
      tolerance_indicator_audit,
      file.path(output_dir, "mc_numerical_tolerance_indicator_audit.csv"),
      row.names = FALSE
    )
    if (any(!tolerance_indicator_audit$pass_zero_for_homothetic, na.rm = TRUE)) {
      warning(
        "A homothetic share indicator remained positive after applying mc_numerical_tol. ",
        "Inspect mc_numerical_tolerance_indicator_audit.csv."
      )
    }
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    mc_plot_data <- subset(
      replication_level_key_summary,
      comparison == "Unit_scale_only_in_Stone"
    )
    if (nrow(mc_plot_data) > 0) {
      scenario_labels <- c(
        normalized_baseline = "Normalized baseline",
        unit_scale_moderate = "Fixed scale: moderate",
        unit_scale_large = "Fixed scale: large",
        unit_scale_mixed = "Fixed scale: mixed signs",
        local_unit_error = "Local-unit error only",
        quality_reporting_only = "Quality/reporting only",
        combined_unit_quality = "Fixed scale plus other contamination",
        high_beta_unit_scale = "Moderate scale, larger expenditure coefficients",
        fixed_scale_only_mixed = "Fixed scale only: mixed signs",
        homothetic_fixed_scale_mixed = "Homothetic null check (beta = 0, sigma_u = 0)"
      )
      mc_plot_data$scenario_label <- unname(
        scenario_labels[as.character(mc_plot_data$scenario)]
      )
      mc_plot_data$scenario_label <- factor(
        mc_plot_data$scenario_label,
        levels = rev(unname(scenario_labels[names(scenarios)]))
      )
      mc_plot_data$elasticity_label <- ifelse(
        mc_plot_data$type == "cross",
        "Cross-price elasticities",
        "Own-price elasticities"
      )
      mc_plot <- ggplot2::ggplot(
        mc_plot_data,
        ggplot2::aes(
          x = scenario_label, y = mean, ymin = p025, ymax = p975
        )
      ) +
        ggplot2::geom_pointrange() +
        ggplot2::coord_flip() +
        ggplot2::facet_wrap(~ elasticity_label, scales = "free_x") +
        ggplot2::labs(
          x = NULL,
          y = paste0(
            "Mean absolute elasticity difference\n",
            "between LA(P, UVraw) and LA(P, UVnorm)"
          )
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(
            margin = ggplot2::margin(t = 6)
          ),
          plot.margin = ggplot2::margin(t = 5.5, r = 8, b = 5.5, l = 8)
        )

      ggplot2::ggsave(
        file.path(output_dir, "fig1.pdf"),
        mc_plot, width = 8.0, height = 5.2, units = "in"
      )
      if (isTRUE(config$save_png)) {
        ggplot2::ggsave(
          file.path(output_dir, "fig1.png"),
          mc_plot, width = 8.0, height = 5.2, units = "in", dpi = 300
        )
      }
    }
  }

  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  message("Monte Carlo completed: ", output_dir)
  invisible(raw_results)
}

# Execution
run_analysis <- function(settings, root = analysis_root) {
  config <- make_config(settings, root)
  check_inputs(config)
  dirs <- c(config$output_dir, config$figure_dir)
  if (config$mode %in% c("empirical", "all")) {
    dirs <- c(dirs, config$input_cache_dir, config$empirical_dir,
              config$robustness_dir, config$audit_dir)
  }
  if (config$mode %in% c("monte_carlo", "all")) dirs <- c(dirs, config$mc_output_dir)
  for (path in dirs) {
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) {
      stop("Cannot create output directory: ", path, call. = FALSE)
    }
  }
  started <- proc.time()[["elapsed"]]
  status <- "incomplete"
  on.exit({
    if (status != "complete") {
      write_run_metadata(config, status, proc.time()[["elapsed"]] - started)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  if (config$mode %in% c("empirical", "all")) {
    message("Building empirical inputs from raw survey files.")
    build_input_data(config)
    prepared <- rebuild_measurement_units(config)
    message("Estimating demand systems and paired bootstrap.")
    empirical <- estimate_empirical(config, prepared)
    if (empirical$bootstrap_design$observed_number_of_strata[1L] != 32L) {
      stop("The replication requires 32 district strata; inspect empirical_bootstrap_design.csv.", call. = FALSE)
    }
    rm(prepared, empirical)
  }
  if (config$mode %in% c("monte_carlo", "all")) {
    message("Running Monte Carlo simulations.")
    run_monte_carlo(config)
    for (extension in c("pdf", "png")) {
      source <- file.path(config$mc_output_dir, paste0("fig1.", extension))
      if (file.exists(source)) {
        file.copy(source, file.path(config$figure_dir, basename(source)), overwrite = FALSE)
      }
    }
  }
  write_run_metadata(config, "complete", proc.time()[["elapsed"]] - started)
  paths <- list.files(config$output_dir, recursive = TRUE, full.names = TRUE)
  paths <- paths[!dir.exists(paths)]
  inventory <- data.frame(
    file = substring(paths, nchar(config$output_dir) + 2L),
    bytes = file.info(paths)$size
  )
  utils::write.csv(inventory, file.path(config$output_dir, "output_manifest.csv"), row.names = FALSE)
  status <- "complete"
  message("Completed. Outputs: ", config$output_dir)
  invisible(config$output_dir)
}

run_analysis(settings)
