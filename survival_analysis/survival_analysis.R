packages <- c("survival", "survminer", "dplyr", "ggplot2", "readr", "readxl", "cowplot")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)}}


# Settings

syng_metadata <- "../data/clinical/SYNG-BC1_curated_export_20260423.xlsx"

X_AXIS_MAX <- 72
BREAK_TIME_BY <- 4
SHOW_CONFIDENCE_INTERVALS <- FALSE
ADMINISTRATIVE_CENSORING <- FALSE

COL_NEG <- "#1f5a91"
COL_POS <- "#F22929"

analyses <- list(
  list(
    name = "SBS3", file = "SBS3/SBS3_patient_level.csv", status_col = "SBS3_status",
    negative_value = "SBS3 negative", positive_value = "SBS3 positive",
    negative_label = "SBS3 negative", positive_label = "SBS3 positive",
    output_dir = "SBS3/SBS3_OS_results", prefix = "SBS3",
    inset_xmin = 42, inset_xmax = 72, inset_ymin = 0.66, inset_ymax = 0.98,
    legend_y = 0.025, legend_x = 0.02
  ),
  list(
    name = "SBS3 + ID6", file = "SBS3_ID6/SBS3_ID6_patient_level.csv", status_col = "HRD_status",
    negative_value = "HRD negative", positive_value = "HRD positive",
    negative_label = "SBS3/ID6 negative", positive_label = "SBS3/ID6 positive",
    output_dir = "SBS3_ID6/SBS3_ID6_OS_results", prefix = "SBS3_ID6",
    inset_xmin = 36, inset_xmax = 74, inset_ymin = 0.70, inset_ymax = 0.99,
    legend_y = 0.025, legend_x = 0.02
  ),
  list(
    name = "SBS6 + SBS15 + SBS44", file = "SBS6_15_44/SBS6_15_44_patient_level.csv", status_col = "MMR_status",
    negative_value = "SBS6/15/44 negative", positive_value = "SBS6/15/44 positive",
    negative_label = "SBS6/15/44 negative", positive_label = "SBS6/15/44 positive",
    output_dir = "SBS6_15_44/SBS6_15_44_OS_results", prefix = "SBS6_SBS15_SBS44",
    inset_xmin = 33, inset_xmax = 74, inset_ymin = 0.70, inset_ymax = 0.995,
    legend_y = 0.025, legend_x = 0.01))


# Survival data

patients <- read_excel(syng_metadata, sheet = "Patients") %>%
  mutate(
    patient_id = trimws(as.character(patient_id)),
    last_FU_OS_days = as.numeric(last_FU_OS_days),
    OS_EVENT = case_when(
      as.character(last_FU_OS_event) %in% c("1", "1.0") ~ 1,
      as.character(last_FU_OS_event) %in% c("0", "0.0") ~ 0,
      tolower(trimws(as.character(last_FU_OS_event))) == "dead" ~ 1,
      tolower(trimws(as.character(last_FU_OS_event))) %in% c("alive", "lost to follow-up", "lost to follow up") ~ 0,
      TRUE ~ NA_real_
    ),
    OS_MONTHS = last_FU_OS_days / 30.4375
  ) %>%
  select(patient_id, last_FU_OS_days, OS_MONTHS, OS_EVENT)

if (anyDuplicated(patients$patient_id)) {
  stop("Patients sheet contains duplicated patient IDs.")}

cat("\nPatients:", nrow(patients))
cat("\nKnown OS duration:", sum(!is.na(patients$OS_MONTHS)))
cat("\nKnown OS event:", sum(!is.na(patients$OS_EVENT)), "\n")


# Helper functions

format_p <- function(p) {
  if (p < 0.0001) return("< 0.0001")
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  sprintf("%.3f", p)}

format_median <- function(row) {
  if (nrow(row) == 0 || is.na(row$median)) return("NR")
  sprintf("%.1f\n(%.1f-%.1f)", row$median, row$lower, row$upper)}


# OS analysis

run_os_analysis <- function(config) {
  
  cat("\n", config$name, "\n")
  
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  
  signature_data <- read_csv(config$file, show_col_types = FALSE) %>%
    mutate(patient_id = trimws(as.character(patient_id))) %>%
    filter(.data[[config$status_col]] %in% c(config$negative_value, config$positive_value))
  
  if ("cohort" %in% names(signature_data) && any(signature_data$cohort != "Synergy")) {
    stop(paste0(config$name, ": patient-level file contains patients outside Synergy."))}
  
  if (anyDuplicated(signature_data$patient_id)) {
    stop(paste0(config$name, ": patient-level file contains duplicated patients."))}
  
  cat("\nGroups before adding survival data:\n")
  print(table(signature_data[[config$status_col]]))
  
  
  # Merge survival data
  
  df <- signature_data %>%
    inner_join(patients, by = "patient_id", relationship = "one-to-one") %>%
    mutate(OS_VALID = !is.na(OS_MONTHS) & OS_MONTHS > 0 & !is.na(OS_EVENT))
  
  cat("\nPatients before OS filtering:", nrow(df))
  cat("\nPatients with valid OS:", sum(df$OS_VALID))
  cat("\nPatients excluded:", sum(!df$OS_VALID), "\n")
  
  if (any(!df$OS_VALID)) {
    print(df %>% filter(!OS_VALID) %>% select(patient_id, last_FU_OS_days, OS_MONTHS, OS_EVENT))}
  
  
  # Final dataset
  
  os_df <- df %>%
    filter(OS_VALID) %>%
    mutate(
      GROUP = factor(.data[[config$status_col]], levels = c(config$negative_value, config$positive_value)),
      KM_GROUP = GROUP,
      original_time = OS_MONTHS,
      original_event = as.integer(OS_EVENT))
  
  cat("\nFinal groups:\n")
  print(table(os_df$GROUP))
  
  cat("\nEvents by group:\n")
  print(table(os_df$GROUP, os_df$original_event))
  
  
  # Administrative censoring
  
  if (ADMINISTRATIVE_CENSORING) {
    os_df <- os_df %>%
      mutate(
        surv_time = pmin(original_time, X_AXIS_MAX),
        surv_event = ifelse(original_time <= X_AXIS_MAX, original_event, 0)
      )
  } else {
    os_df <- os_df %>%
      mutate(surv_time = original_time, surv_event = original_event)}
  
  
  # Kaplan-Meier, log-rank and Cox
  
  fit <- survfit(Surv(surv_time, surv_event) ~ KM_GROUP, data = os_df)
  logrank <- survdiff(Surv(surv_time, surv_event) ~ KM_GROUP, data = os_df)
  p_value <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
  
  cox_model <- coxph(Surv(surv_time, surv_event) ~ GROUP, data = os_df)
  cox_summary <- summary(cox_model)
  ph_test <- cox.zph(cox_model)
  median_os <- surv_median(fit)
  
  cat("\nCox model:\n")
  print(cox_summary)
  
  cat("\nProportional hazards assumption:\n")
  print(ph_test)
  
  cat("\nMedian OS:\n")
  print(median_os)
  
  
  # Plot values
  
  n_neg <- sum(os_df$GROUP == config$negative_value)
  n_pos <- sum(os_df$GROUP == config$positive_value)
  
  med_neg <- median_os %>% filter(strata == paste0("KM_GROUP=", config$negative_value))
  med_pos <- median_os %>% filter(strata == paste0("KM_GROUP=", config$positive_value))
  
  med_neg_label <- format_median(med_neg)
  med_pos_label <- format_median(med_pos)
  
  hr <- cox_summary$coefficients[1, "exp(coef)"]
  hr_low <- cox_summary$conf.int[1, "lower .95"]
  hr_up <- cox_summary$conf.int[1, "upper .95"]
  cox_p <- cox_summary$coefficients[1, "Pr(>|z|)"]
  
  hr_label <- sprintf("%.2f (%.2f-%.2f)", hr, hr_low, hr_up)
  logrank_label <- format_p(p_value)
  
  
  # Survival at 24 and 36 months
  
  landmark_times <- c(24, 36)
  landmark_sum <- summary(fit, times = landmark_times, extend = TRUE)
  
  landmark_df <- data.frame(
    strata = sub("^KM_GROUP=", "", landmark_sum$strata),
    time = landmark_sum$time,
    surv = landmark_sum$surv)
  
  get_surv <- function(group, time_point) {
    value <- landmark_df$surv[landmark_df$strata == group & landmark_df$time == time_point]
    if (length(value) == 0) return(NA_real_)
    value[1]}
  
  s24_neg <- get_surv(config$negative_value, 24)
  s24_pos <- get_surv(config$positive_value, 24)
  s36_neg <- get_surv(config$negative_value, 36)
  s36_pos <- get_surv(config$positive_value, 36)
  
  
  # Kaplan-Meier plot
  
  km <- ggsurvplot(
    fit, data = os_df,
    conf.int = SHOW_CONFIDENCE_INTERVALS,
    risk.table = TRUE,
    risk.table.height = 0.1,
    censor = TRUE,
    censor.shape = 3,
    censor.size = 2.2,
    palette = c(COL_NEG, COL_POS),
    legend.title = "",
    legend.labs = c(
      paste0(config$negative_label, " (n = ", n_neg, ")"),
      paste0(config$positive_label, " (n = ", n_pos, ")")
    ),
    legend = "bottom",
    xlim = c(0, X_AXIS_MAX),
    ylim = c(0, 1),
    break.time.by = BREAK_TIME_BY,
    xlab = "Time (months)",
    ylab = "OS probability (%)",
    pval = FALSE,
    risk.table.y.text = TRUE,
    risk.table.y.text.col = FALSE,
    risk.table.title = "Patients still at risk:",
    risk.table.fontsize = 5,
    ggtheme = theme_classic(base_size = 14))
  
  main_plot <- km$plot +
    scale_y_continuous(
      breaks = seq(0, 1, 0.1),
      labels = function(x) paste0(x * 100, "%"),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      breaks = seq(0, X_AXIS_MAX, BREAK_TIME_BY),
      limits = c(0, X_AXIS_MAX),
      expand = c(0.01, 0))
  
  
  # Landmark lines
  
  if (!is.na(s24_neg) && !is.na(s24_pos)) {
    main_plot <- main_plot +
      geom_segment(
        aes(x = 24, xend = 24, y = 0, yend = max(s24_neg, s24_pos)),
        linetype = "dashed", linewidth = 0.5, colour = "grey55"
      ) +
      geom_segment(
        aes(x = 0, xend = 24, y = s24_neg, yend = s24_neg),
        linetype = "dashed", linewidth = 0.5, colour = "grey65"
      ) +
      geom_segment(
        aes(x = 0, xend = 24, y = s24_pos, yend = s24_pos),
        linetype = "dashed", linewidth = 0.5, colour = "grey65")}
  
  if (!is.na(s36_neg) && !is.na(s36_pos)) {
    main_plot <- main_plot +
      geom_segment(
        aes(x = 36, xend = 36, y = 0, yend = max(s36_neg, s36_pos)),
        linetype = "dashed", linewidth = 0.5, colour = "grey55"
      ) +
      geom_segment(
        aes(x = 0, xend = 36, y = s36_neg, yend = s36_neg),
        linetype = "dashed", linewidth = 0.5, colour = "grey65"
      ) +
      geom_segment(
        aes(x = 0, xend = 36, y = s36_pos, yend = s36_pos),
        linetype = "dashed", linewidth = 0.5, colour = "grey65")}
  
  
  # Censoring legend
  
  censor_legend_data <- data.frame(x = 0, y = 0, legend_type = "Censored")
  
  main_plot <- main_plot +
    geom_point(
      data = censor_legend_data, aes(x = x, y = y, shape = legend_type),
      inherit.aes = FALSE, colour = "black", alpha = 0, show.legend = TRUE
    ) +
    scale_shape_manual(name = NULL, values = c("Censored" = 3)) +
    guides(
      colour = guide_legend(order = 2, override.aes = list(shape = NA)),
      shape = guide_legend(order = 1, override.aes = list(alpha = 1, size = 4, colour = "black"))
    ) +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_text(size = 14, margin = margin(t = 12)),
      axis.title.y = element_text(size = 14, margin = margin(t = 12)),
      axis.text = element_text(size = 12),
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.5),
      legend.position = c(config$legend_x, config$legend_y),
      legend.justification = c(0, 0),
      legend.box = "vertical",
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.spacing.y = grid::unit(-0.55, "cm"),
      legend.key.width = grid::unit(1.2, "cm"),
      legend.text = element_text(size = 14),
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5))
  
  
  # Risk table
  
  risk_table <- km$table +
    scale_x_continuous(
      breaks = seq(0, X_AXIS_MAX, BREAK_TIME_BY),
      limits = c(0, X_AXIS_MAX),
      expand = c(0.01, 0)
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 11, face = "plain", hjust = 0),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12, margin = margin(r = 10)),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 5, 2, 5))
  
  
  # Statistics inset
  
  inset_table <- ggplot() +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    annotate("segment", x = 0.001, xend = 0.98, y = 0.75, yend = 0.75, linewidth = 0.5) +
    annotate("segment", x = 0.001, xend = 0.98, y = 0.45, yend = 0.45, linewidth = 0.5) +
    annotate("segment", x = 0.001, xend = 0.98, y = 0.25, yend = 0.25, linewidth = 0.5) +
    annotate("text", x = 0.56, y = 0.93, label = paste0(config$negative_label, "\n(n = ", n_neg, ")"), size = 4) +
    annotate("text", x = 0.86, y = 0.93, label = paste0(config$positive_label, "\n(n = ", n_pos, ")"), size = 4) +
    annotate("text", x = 0.01, y = 0.60, hjust = 0, label = "Median (95% CI),\nmonths", size = 4) +
    annotate("text", x = 0.56, y = 0.60, label = med_neg_label, size = 4) +
    annotate("text", x = 0.86, y = 0.60, label = med_pos_label, size = 4) +
    annotate("text", x = 0.01, y = 0.35, hjust = 0, label = "Hazard ratio (95% CI)", size = 4) +
    annotate("text", x = 0.70, y = 0.35, label = hr_label, size = 4) +
    annotate("text", x = 0.01, y = 0.16, hjust = 0, label = "Log-rank p-value", size = 4) +
    annotate("text", x = 0.70, y = 0.16, label = logrank_label, size = 4)
  
  main_plot <- main_plot +
    annotation_custom(
      grob = ggplotGrob(inset_table),
      xmin = config$inset_xmin, xmax = config$inset_xmax,
      ymin = config$inset_ymin, ymax = config$inset_ymax)
  
  final_plot <- cowplot::plot_grid(
    main_plot, risk_table,
    ncol = 1, rel_heights = c(4.3, 0.85), align = "v", axis = "lr")
  
  print(final_plot)
  
  
  # Save plot
  
  ggsave(
    file.path(config$output_dir, paste0("KM_OS_", config$prefix, ".png")),
    plot = final_plot, width = 12.5, height = 8.2, dpi = 300, bg = "white")
  
  
  # Save results
  
  group_summary <- os_df %>%
    group_by(GROUP) %>%
    summarise(n = n(), deaths = sum(surv_event == 1), censored = sum(surv_event == 0), .groups = "drop")
  
  logrank_result <- data.frame(
    comparison = paste(config$negative_label, "vs", config$positive_label),
    p_value = p_value, n_total = nrow(os_df), n_negative = n_neg, n_positive = n_pos)
  
  cox_result <- data.frame(
    comparison = paste(config$positive_label, "vs", config$negative_label),
    hazard_ratio = hr, CI_lower_95 = hr_low, CI_upper_95 = hr_up, p_value = cox_p)
  
  ph_result <- as.data.frame(ph_test$table)
  ph_result$term <- rownames(ph_result)
  rownames(ph_result) <- NULL
  
  write.csv(os_df, file.path(config$output_dir, paste0(config$prefix, "_OS_final_dataset.csv")), row.names = FALSE)
  write.csv(group_summary, file.path(config$output_dir, paste0(config$prefix, "_OS_group_summary.csv")), row.names = FALSE)
  write.csv(logrank_result, file.path(config$output_dir, paste0(config$prefix, "_OS_logrank_result.csv")), row.names = FALSE)
  write.csv(cox_result, file.path(config$output_dir, paste0(config$prefix, "_OS_cox_result.csv")), row.names = FALSE)
  write.csv(median_os, file.path(config$output_dir, paste0(config$prefix, "_OS_KM_median.csv")), row.names = FALSE)
  write.csv(ph_result, file.path(config$output_dir, paste0(config$prefix, "_OS_PH_test.csv")), row.names = FALSE)
  
  cat("\nGroup summary:\n")
  print(group_summary)
  cat("\nLog-rank p-value:", p_value, "\n")
  cat("\nCox result:\n")
  print(cox_result)}


# Run analyses

for (config in analyses) {
  run_os_analysis(config)}