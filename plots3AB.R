library(tidyverse)
library(purrr)
library(patchwork)

files <- c(
  full       = "test_r2_existing_models_long.csv",
  no_nkp     = "test_r2_no_nkp_long.csv",
  no_nkp_56  = "test_r2_no_nkp_nkcd56_long.csv",
  no_all_nk  = "test_r2_no_all_nk_long.csv"
)

dat <- imap(files, ~ read_csv(.x, show_col_types = FALSE) %>%
              mutate(run = .y)) %>%
  list_rbind()

dat$run <- factor(dat$run, levels = names(files))
NK <- c("NK", "NK_CD56bright", "NK_Proliferating")

run_labels_3 <- c(
  full       = "All contexts",
  no_nkp     = "−NK_Prolif",
  no_nkp_56  = "−NK_Prolif,\n−NK_CD56br",
  no_all_nk  = "−all NK"
)

panel_a <- dat %>%
  filter(method == "Full_Model", r2 > 0) %>%
  mutate(group = if_else(context %in% NK, "NK contexts", "Other contex")) %>%
  count(run, group) %>%
  ggplot(aes(run, n, fill = group)) +
  geom_col(width = 0.7) +
  scale_x_discrete(labels = run_labels_3, drop = FALSE) +
  scale_fill_manual(values = c("NK contexts" = "#d55e00",
                               "Other contex" = "#999999")) +
  labs(x = "Ablation run",
       y = "Predictions with test R-squared > 0",
       fill = NULL) +
  theme_minimal() +
  theme(plot.tag = element_text(size = 16, face = "bold"))

panel_b <- dat %>%
  filter(context %in% NK, method == "Full_Model") %>%
  group_by(run, context) %>%
  summarise(med_r2 = median(r2), .groups = "drop") %>%
  ggplot(aes(run, med_r2, color = context, group = context)) +
  geom_line() +
  geom_point(size = 2.5) +
  scale_x_discrete(labels = run_labels_3, drop = FALSE) +
  coord_cartesian(ylim = c(0.10, 0.18)) +
  labs(x = "Ablation run",
       y = "Median test R-squared (Full Model)",
       color = "NK context") +
  theme_minimal() +
  theme(plot.tag = element_text(size = 16, face = "bold"))

fig3_combined <- (panel_a + labs(tag = "A")) | (panel_b + labs(tag = "B"))

fig3_combined
ggsave("fig3_combined.png", fig3_combined, width = 12, height = 4.5, dpi = 300)