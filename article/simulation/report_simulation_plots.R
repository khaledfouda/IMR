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
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = seq(.02, 0.3, 0.04)) +
  facet_wrap(~ metric_lab, labeller = label_parsed, ncol=3, scales="free_y") +
  labs(x = "Observed Rate (Training)", y = NULL, color = "Model", fill = "Model") +
  theme_minimal(base_size = 10) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    strip.text = element_text(face = "bold")
  ) -> sim2.g; sim2.g

ggsave("./article/simulation/data/sim2_plot1.png", sim2.g, width = 320/25.4, height = 150/25.4, dpi = 600)
#--------------------------------------------------------------------------------




