library(dplyr)


result_df = replicate_25_2$summary


# numeric summary
summary_df <- result_df %>%
  group_by(term, scenario_id, scenario_label, n_total, n_sites) %>%
  summarise(
    agreement_rate = mean(agreement),
    sd_agreement = sd(agreement),
    mean_diff = mean(diff, na.rm = TRUE),
    sd_diff = sd(diff, na.rm = TRUE),
    mean_log_diff = mean(log_diff, na.rm = TRUE),
    sd_log_diff = sd(log_diff, na.rm = TRUE),
    mean_p = mean(p, na.rm = TRUE),
    mean_p_meta = mean(p_meta, na.rm = TRUE),
    .groups = "drop"
  )

