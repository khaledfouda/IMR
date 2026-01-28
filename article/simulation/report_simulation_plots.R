library(devtools)
#clean_dll(); Rcpp::compileAttributes(); document()
load_all()
# devtools::uninstall(); devtools::install()
require(tidyverse)
require(magrittr)
source("./article/simulation/generate_simu_dat.R")
source("./other_models/SoftImpute_cv.R")
source("./other_models/MCCI.R")


results1 <- readRDS("article/simulation/data/sim1_res.rds") %>%
  mutate(metric = if_else(test==test_rel, "Rel.RMSE", "RMSE"))
results2 <- readRDS("article/simulation/data/sim2_res.rds")

results1 %>%
  filter(metric == "RMSE") %>%
  group_by(dim, model, metric) %>%
  summarise_all(mean)

#-- Simulation 1 - table 1 -----------------------
compute.mean.sd <- function(x){
  s <- sd(x)
  if(is.na(s)) "" else paste0(round(mean(x),3)," (", round(s,3) ,")")
}

results1 %>%
  filter(metric != "RMSE") %>%
  dplyr::select(-gamma, -metric, -test_rel, -theta) %>%
  arrange(dim, model) %>%
  dplyr::group_by(dim, model) %>%
  summarize_all(compute.mean.sd)  %>%
  pivot_longer(c("beta", "M",  "train", "test", "rank"),
               names_to = c("metric")) %>%
  mutate(value = if_else(model == "SI" & metric == "M","",value)) %>%
  mutate(model = if_else(model == "SI", "SoftImpute", model)) ->
  sim1.tab


require(kableExtra)

metric_labels <- c(
  beta = "RRMSE($\\beta$)",
  M    = "RRMSE($M$)",
  #theta    = "RMSE($\\Theta$)",
  train = "RRMSE(train)",
  test = "RRMSE(test)",
  rank       = "Rank"
)

# (Optional) method display order (put yours here to match the paper)
method_order <- c(
  "IMR","SoftImpute",  "MCCI"
)

wide <- sim1.tab %>%
  mutate(metric = recode(metric, !!!metric_labels)) %>%
  mutate(model = factor(model, levels = method_order)) %>%
  arrange(dim, model) %>%
  select(dim, model, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  as.data.frame()

# How many rows per panel (for pack_rows)
panel_counts <- wide %>% count(dim, name = "rows")
panel_index  <- stats::setNames(panel_counts$rows,
                                paste0("n = m = ", panel_counts$dim))
rows_to_bold <- c(1, 4, 7, 10)
bold_vector <- seq_len(nrow(wide)) %in% rows_to_bold
# Build table
kbl(
  wide %>% dplyr::select(-dim),
  format    = "html",
  booktabs  = TRUE,
  escape    = FALSE,
  label =  "tab:sim1",
  position  = "htbp",
  # col.names = c("Method", setdiff(names(wide), c("n", "method"))),
  align     = c("l", rep("r", ncol(wide) - 2)),
  caption   = paste(
    "Empirical relative root mean square errors (RRMSEs),",
    "estimated ranks, and their standard deviations (in parentheses)",
    "under model $\\Theta=X\\beta+M$ and with dimensions $(n,m)=(400,400),(600,600),(800,800),(1000,1000)$,",
    "$80\\%$ missingness rate, and the true rank is 14.")
) |>
  pack_rows(index = panel_index, bold = FALSE) |>
  kable_styling(latex_options = c("hold_position", "striped"), font_size = 12) |>
  column_spec(5, color = "black",bold = bold_vector) |>
  column_spec(1, color = "black",bold = TRUE ) |>
  row_spec(1:nrow(wide), color = "black") -> sim1.tbl;sim1.tbl

#==============================================================================

#--- simulation 2 - Plot 1 >>
results2 <- readRDS("article/simulation/data/sim2_res.rds")

results2 %>%
  filter(model == "SI") %>%
  mutate(miss_pct = round(miss_pct, 2)) %>%
  select(rank, miss_pct) %>%
  group_by(miss_pct) %>%
  summarise_all(mean)

results2 %>%
  filter(metric != "RMSE") %>%
  #dplyr::mutate(rank = abs(rank - 15)) %>%
  dplyr::mutate(sparsity = round(miss_pct, 2)) %>%
  dplyr::select(model, sparsity, beta, gamma, M, theta, test, rank) %>%
  dplyr::group_by(model, sparsity) %>%
  dplyr::mutate(M = if_else(model == "SI", NA, M)) %>%
  dplyr::mutate(model = if_else(model == "SI", "SoftImpute", model)) %>%
  dplyr::summarize_all(c(error_mean=mean,error_sd= sd)) %>%
  dplyr::ungroup() %>%
  pivot_longer(-c(model, sparsity),
               names_to = c("metric", "stat"),
               names_pattern = "^(.*)_error_(mean|sd)$",
               values_to = "val") %>%
  pivot_wider(names_from = stat, values_from = val) %>%
  mutate(sparsity = 1-sparsity) %>%
  arrange(model, sparsity) %>%
  mutate(ymin = pmax(0, mean-sd), ymax = mean+sd) ->
  sim2.long


metric_labels <- c(
  beta = "RRMSE(beta)",
  gamma = "RRMSE(Gamma)",
  M    = "RRMSE(M)",
  theta    = "RRMSE(Theta)",
  test = "RRMSE(test)",
  rank       = "Estimated~Rank"
)
sim2.long %<>%
  mutate(metric_lab = factor(metric,
                             levels=names(metric_labels),
                             labels = unname(metric_labels)))
okabe_ito <- c("#56B4E9","#E69F00")
metrics <- unique(sim2.long$metric)
require(scales)
ggplot(sim2.long, aes(x = sparsity, y = mean, color = model, fill = model, group = model)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA) +
  geom_line(size = 1) +
  geom_point(size = 1.6) +
  scale_color_manual(values = okabe_ito) +
  scale_fill_manual(values  = okabe_ito) +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = seq(.05, 0.3, length=6)) +
  facet_wrap(~ metric_lab, labeller = label_parsed, ncol=3, scales="free_y") +
  labs(x = "Observed Rate (Training Size)", y = NULL, color = "Model", fill = "Model") +
  theme_minimal(base_size = 10) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    strip.text = element_text(face = "bold")
  ) -> sim2.g; sim2.g

ggsave("./article/simulation/data/sim2_plot1.png", sim2.g, width = 320/25.4, height = 150/25.4, dpi = 600)
#--------------------------------------------------------------------------------
results3 <- readRDS("article/simulation/data/sim2_2_res.rds")

results3 %>%
  filter(metric == "Rel.RMSE") %>%
  mutate(obs_rate = 1 - round(miss_pct, 2)) %>%
  group_by(model, obs_rate) %>%
  summarise(
    mean_time = mean(time, na.rm = TRUE),
    mean_error = mean(test, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = model,
    values_from = c(mean_time, mean_error)
  ) %>%
  mutate(
    time_ratio = mean_time_IMR / mean_time_SI,
    error_improve = (mean_error_SI - mean_error_IMR) / mean_error_SI
  ) %>%
  select(obs_rate, time_ratio, error_improve) %>%
  pivot_longer(
    cols = c(time_ratio, error_improve),
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(
    measure_label = case_when(
      measure == "time_ratio" ~ "Computational Cost (Times Slower than SoftImpute)",
      measure == "error_improve" ~ "Performance Gain (% Improvement over SoftImpute)"
    )
  ) -> plot_data

ggplot(plot_data, aes(x = obs_rate, y = value)) +
  geom_hline(data = filter(plot_data, measure == "time_ratio"),
             aes(yintercept = 1), linetype = "dashed", color = "grey50") +
  geom_hline(data = filter(plot_data, measure == "error_improve"),
             aes(yintercept = 0), linetype = "dashed", color = "grey50") +
  geom_line(size = 1, color = "#E69F00") + # Okabe-Ito Orange
  geom_point(size = 2, color = "#E69F00") +
  facet_wrap(~ measure_label, scales = "free_y", ncol = 1) +
  scale_x_continuous("Observation Rate (Training Size)", labels = percent_format(accuracy = 1),
                     breaks = seq(0.05, 0.3, length.out=6)) +
  scale_y_continuous(
    name = NULL,
    labels = function(x) {
      ifelse(x < 1 & x > -1, percent(x, accuracy = 1), number(x, accuracy = 0.1, suffix = "x"))
    }
  ) +

  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95")
  ) -> diff_plot; diff_plot

ggsave("./article/simulation/data/sim2_plot2.png", diff_plot, width = 320/25.4, height = 150/25.4, dpi = 600)


