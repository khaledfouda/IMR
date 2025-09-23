library(devtools)
load_all()

require(tidyverse)
library(magrittr)
source("./notes/generate_simu_dat.R")
source("./other_models/SoftImpute_cv.R")
source("./other_models/MCCI.R")
source("./other_models/Naive.R")

sim1_res <- function(dat, fit, name="", ortho=TRUE, error_metric=IMR::error_metric$rmse){
  # prepare data : we need values for: M, beta, theta
  # expect fit$ to contain (u, d, and v) or (M)

  has.beta <- "beta" %in% names(fit) & ! is.null(fit$beta)
  has.gamma <- "gamma" %in% names(fit) & ! is.null(fit$gamma)
  has.M    <- "M" %in% names(fit) | all(c("u","d","v") %in% names(fit))
  has.intercept <- "beta0" %in% names(fit) & !is.null(fit$beta0)
  estimates <- 0
  out <- data.frame(model=name, beta=NA, gamma=NA, M=NA, theta=NA, test=NA, rank=NA)
  # check M
  if(all(c("u", "d", "v") %in% names(fit)))
    fit$M <- fit$u %*% (fit$d * t(fit$v))
  if(has.M) {
    estimates <- fit$M
    out$M <- error_metric(fit$M, dat$M)
  }
  if(has.beta){
    fit$beta <- solve(dat$fit_data$X$R) %*% fit$beta
    estimates <- estimates + dat$X %*% fit$beta
    out$beta <- error_metric(fit$beta, dat$beta)
  }
  if(has.gamma){
    fit$gamma <- fit$gamma %*% solve(t(dat$fit_data$Z$R))
    estimates <- estimates + fit$gamma %*% t(dat$Z)
    out$gamma <- error_metric(fit$gamma, dat$gamma)
  }
  if(has.intercept){
    estimates <- estimates + fit$beta0 %*% matrix(1,1,ncol(dat$Y))
  }


  stopifnot(all(estimates!=0))
  out$theta <- error_metric(estimates, dat$theta)
  test.obs <- dat$Y == 0
  out$test <- IMR::error_metric$rel.rmse(estimates[test.obs], dat$theta[test.obs])
  out$rank <- qr(estimates)$rank
  return(out)
}


one_loop_fit <- function(dat, hpar = IMR::get_imr_default_hparams(), seed=2025){


 # fit.imr <- IMR::imr.cv(dat$fit_data$Y_full, X=dat$fit_data$X$Q,
 #                      Z = dat$fit_data$Z$Q,intercept_row = F,
 #                      hpar = hpar, seed = seed, ls_initial = FALSE,
 #                      intercept_col = F, verbose=1)

  fit.nlrr <-  IMR:::nlrr.cv(dat$fit_data$Y_full, X=dat$fit_data$X$Q,
                           Z = dat$fit_data$Z$Q,intercept_row = F,
                           hpar = hpar, seed = seed,
                           intercept_col = F, verbose=1)

  fitsi <- simpute.cv(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$Y_full,
                      trace=T,
                      tol = hpar$M$early.stopping, seed = seed,
                      n.lambda = hpar$M$n.lambda
                      #, n.lambda=80, rank.init = 2,
                      #rank.step = 1,
                      #lambda0_fun = IMR::get_lambda_M_max
  )
  fitmmci <- MCCI.cv(dat$Y, dat$X, dat$mask, numCores = 9, seed = seed)
  fitnai <- naive.fit(dat$Y, dat$X)

  rbind(
  #sim1_res(dat, fit.imr$fit, "IMR"),
  sim1_res(dat, fit.nlrr$fit, "IMR"),
  sim1_res(dat, fitsi$fit, "SoftImpute"),
  sim1_res(dat, fitmmci$fit, "MCCI"),
  sim1_res(dat, fitnai, "Naive")
  )
}





future::plan(future::sequential)
future::plan(future::multisession, workers = 9)

hpar <- IMR::get_imr_default_hparams()

# hpar$M$lambda_factor <- 1
hpar$M$n.lambda <- 30
# hpar$M$rank.init <- 2
# hpar$M$rank.step <- 2
# hpar$M$rank.max  <- 15
# hpar$M$lambda_max <- 80
# hpar$beta$n.lambda <- 10
# hpar$gamma$n.lambda <- 10
hpar$M$early.stopping <- 3

results5 <- data.frame()

for(i in 7:30){
  seed = 2025 + i
  dat <-
    generate_simulated_data(1000, 1000, 10, 20, 0, 0.8,
                            sparsity_beta = .5, sparsity_gamma = 0.0,
                            prepare_for_fitting = T,mv_coeffs = T,seed = seed)
  start <- Sys.time()
  results5 %<>% rbind(
  one_loop_fit(dat, hpar, seed)
  )
  results5 %>%
    dplyr::group_by(model) %>%
    dplyr::summarise_all(mean) %>%
    print()
  print(Sys.time() - start)
  print(i)

}

saveRDS(results5, "./notes/saved/results51000.rds")

results5 %>%
  dplyr::group_by(model) %>%
  summarize_all(compute.mean.sd )

#----- combine results
compute.mean.sd <- function(x){
  s <- sd(x)
  if(is.na(s)) "" else paste0(round(mean(x),3)," (", round(s,3) ,")")
}

results3 <- readRDS("./notes/saved/results2800.rds")
results2 <- readRDS("./notes/saved/results2600.rds")
results4 <- readRDS("./notes/saved/results4400.rds")
results5 <- readRDS("./notes/saved/results51000.rds")

results3 %>%
  dplyr::select(-gamma) %>%
  dplyr::group_by(model) %>%
  summarize_all(compute.mean.sd ) %>%
  mutate(n = 800) %>%
  pivot_longer(c("beta", "M", "theta", "test", "rank"),
               names_to = c("metric")) %>%
  rbind(
       results2 %>%
         dplyr::select(-gamma) %>%
         dplyr::group_by(model) %>%
         summarize_all(compute.mean.sd ) %>%
         mutate(n = 600) %>%
         pivot_longer(c("beta", "M", "theta", "test", "rank"), names_to = c("metric"))
       ) %>%
  rbind(
    results4 %>%
      dplyr::select(-gamma) %>%
      dplyr::group_by(model) %>%
      summarize_all(compute.mean.sd ) %>%
      mutate(n = 400) %>%
      pivot_longer(c("beta", "M", "theta", "test", "rank"), names_to = c("metric"))
  ) %>%
  rbind(
    results5 %>%
      dplyr::select(-gamma) %>%
      dplyr::group_by(model) %>%
      summarize_all(compute.mean.sd ) %>%
      mutate(n = 1000) %>%
      pivot_longer(c("beta", "M", "theta", "test", "rank"), names_to = c("metric"))
  )-> sim1.tab

# saveRDS(tab, "./notes/saved/tab_sim1.rds")
sim1.tab <- readRDS("./notes/saved/tab_sim1.rds")
require(kableExtra)



metric_labels <- c(
  beta = "RMSE($\\beta$)",
  M    = "RMSE($M$)",
  theta    = "RMSE($\\Theta$)",
  test = "Test error",
  rank       = "Rank($\\Theta$)"
)

# (Optional) method display order (put yours here to match the paper)
method_order <- c(
  "IMR", "SoftImpute", "MCCI",  "Naive"
)

wide <- sim1.tab %>%
  mutate(metric = recode(metric, !!!metric_labels)) %>%
  mutate(model = factor(model, levels = method_order)) %>%
  arrange(n, model) %>%
  select(n, model, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)

# How many rows per panel (for pack_rows)
panel_counts <- wide %>% count(n, name = "rows")
panel_index  <- stats::setNames(panel_counts$rows,
                                paste0("$n = m = ", panel_counts$n, "$"))

# Build table
kbl(
  wide %>% select(-n),
  format    = "latex",
  booktabs  = TRUE,
  escape    = FALSE,  # keep LaTeX math in headers
 # col.names = c("Method", setdiff(names(wide), c("n", "method"))),
  align     = c("l", rep("r", ncol(wide) - 2)),
  caption   = paste(
                    "Empirical root mean square errors (RMSEs), test error,",
                    "estimated ranks, and their standard errors (in parenthes)",
  "under model $\\Theta=X\\beta+M$ and with $(n,m)=(400,400),(600,600),(800,800),(1000,1000)$,",
  "rank$(\\Theta)=30$, $r = 10$, $p= 20$, 20% observation rate of $\\Theta$,",
  "and 50% of $\\beta$ is 0, at random.")
) |>
  pack_rows(index = panel_index) |>
  kable_styling(latex_options = c("hold_position", "striped"), font_size = 12) |>
  row_spec(0, color = "black") |>
  row_spec(1:nrow(wide), color = "black") -> sim1.tbl;sim1.tbl


osr <- function(n,m,r, pi=NULL, rho=NULL){
  if(is.null(pi)){
    stopifnot(!is.null(rho))
    pi = (rho * r * (n+m-r))/(n*m)
    return(pi)
  }else if(is.null(rho)){
    stopifnot(! is.null(pi))
    rho = (pi * (n*m)) / (r*(n+m-r))
    return(rho)
  }
  stop("something is wrong")
}

for(n in seq(200, 1200,200))
  print(osr(n, n, 10, rho=5))

#==================================================================
#' we now work on second part of the simulation
#' Results needed: Rmse > beta, gamma, M, test
#'                 obtain an average hyperparameters (r, 3 x lambda)
#' Variables: Fix n,m,p,q; change observation rate
#' Models : IMC, SoftImpute
#' Second step: use the fit function only to get speed vs model and observation rate
increase_sparsity <- function(dat, step=0.05){
  current_length <- length(dat$fit_data$Y_full@x)
  current_sparsity <- 1- current_length / length(dat$Y)
  target_sparsity <- step + current_sparsity
  print(paste("Target sparsity is ", target_sparsity))
  stopifnot(target_sparsity < 1)
  extra_nonzero_frac <- step / (1 - current_sparsity)

  to_zero_ind <- sample(1:current_length, extra_nonzero_frac*current_length,replace = F)
  dat$fit_data$Y_full@x[to_zero_ind] <- NA
  dat$fit_data$Y_full %<>% IMR::as.Incomplete()
  #-- we now recreate the train/test splits
  dat$mask <- as.matrix(dat$fit_data$Y_full != 0)
  valid_mask <- MC_train_test_split(dat$mask, 0.2)
  dat$fit_data$train <- dat$fit_data$Y_full * valid_mask
  dat$fit_data$valid <- dat$fit_data$Y_full * (1-valid_mask)
  dat$sparsity <- target_sparsity
  return(dat)
}


lambdas <- data.frame()
sim2.res <- data.frame()
future::plan(future::sequential)
future::plan(future::multisession, workers = 9)

hpar <- IMR::get_imr_default_hparams()
hpar$M$n.lambda <- 40
hpar$M$early.stopping <- 3

for(b in 1:30){
seed = 2025 + b
dat <-
  generate_simulated_data(1000, 1000, 10, 20, 20, 0.7-.02,
                          sparsity_beta = .5, sparsity_gamma = 0.5,
                          prepare_for_fitting = T,mv_coeffs = T,seed = seed)


for(s in 1:15){

dat <- increase_sparsity(dat, .02)


fit.nlrr <-  IMR:::nlrr.cv(dat$fit_data$Y_full, X=dat$fit_data$X$Q,
                           Z = dat$fit_data$Z$Q,intercept_row = F,
                           hpar = hpar, seed = seed,
                           intercept_col = F, verbose=1)

lambdas %<>% rbind(
  data.frame(sparsity=dat$sparsity, n=1000,beta=fit.nlrr$lambda_beta,
             gamma = fit.nlrr$lambda_gamma,
             M = fit.nlrr$lambda_M,
             r = fit.nlrr$rank_M)
)

fitsi <- simpute.cv(dat$fit_data$train, dat$fit_data$valid, dat$fit_data$Y_full,
                    trace=F,
                    tol = hpar$M$early.stopping, seed = seed,
                    n.lambda = hpar$M$n.lambda)

sim2.res %<>% rbind(
  sim1_res(dat, fit.nlrr$fit, "IMR") %>% cbind(sparsity = dat$sparsity, n=1000),
  sim1_res(dat, fitsi$fit, "SoftImpute") %>% cbind(sparsity = dat$sparsity, n=1000)
 )
print(sim2.res)
}
saveRDS(lambdas, "./notes/saved/sim2_parameters.rds")
saveRDS(sim2.res, "./notes/saved/sim2_results.rds")
}
#--------------------------------


lambdas %>%
  mutate(sparsity = round(sparsity, 2)) %>%
  group_by(sparsity) %>%
  summarize_all(mean) %>%
  # summarise_all(function(x) paste(round(mean(x),5),"-",round(sd(x),5))) %>%
  ungroup()



sim2.res %>%
  dplyr::mutate(rank = abs(rank - 50)) %>%
  dplyr::mutate(sparsity = round(sparsity, 2)) %>%
  dplyr::group_by(model, sparsity, n) %>%
  dplyr::summarize_all(c(error_mean=mean,error_sd= sd)) %>%
  dplyr::ungroup() %>%
  select(-n) %>%
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
  beta = "RMSE(beta)",
  gamma = "RMSE(Gamma)",
  M    = "RMSE(M)",
  theta    = "RMSE(Theta)",
  test = "Test~error",
  rank       = "Rank~error"
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
  # facet_wrap(~ factor(metric, levels = metrics, labels = metric_labels[metrics]),
  #            scales = "free_y", ncol = 3) +
  labs(x = "Observation rate", y = "Error", color = "Model", fill = "Model") +
   theme_minimal(base_size = 10) +
  #cowplot::theme_cowplot()+
  theme_bw() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    strip.text = element_text(face = "bold")
  ) -> sim2.g; sim2.g
ggsave("./notes/saved/sim2_plot.png", sim2.g, width = 320/25.4, height = 150/25.4, dpi = 600)

#==============
#-- now the second part of simulation 2: < Time vs fit >


lambdas %>%
  mutate(sparsity = round(sparsity, 2)) %>%
  group_by(sparsity) %>%
  summarize_all(mean) %>%
  # summarise_all(function(x) paste(round(mean(x),5),"-",round(sd(x),5))) %>%
  ungroup() %>%
  select(-n) -> lambdas


future::plan(future::sequential)
future::plan(future::multisession, workers = 9)

time_results <- data.frame()

for(b in 1:30){
  seed = 2025 + b
  dat <-
    generate_simulated_data(1000, 1000, 10, 20, 20, 0.7-.02,
                            sparsity_beta = .5, sparsity_gamma = 0.5,
                            prepare_for_fitting = T,mv_coeffs = T,seed = seed)


  for(s in 1:15){

    dat <- increase_sparsity(dat, .02)
    dat$sparsity <- round(dat$sparsity,2)
    stopifnot(dat$sparsity == lambdas$sparsity[s])

    bench::bench_time(
      IMR::imr.fit(
      Y = dat$fit_data$Y_full,
      X = dat$fit_data$X$Q,
      Z = dat$fit_data$Z$Q,
      r = lambdas$r[s],
      lambda_M = lambdas$M[s],
      lambda_beta = lambdas$beta[s],
      lambda_gamma = lambdas$gamma[s],
      trace = FALSE,
      ls_initial = T
    )) -> imrt


    bench::bench_time(softImpute::softImpute(dat$Y,
                                    rank.max = lambdas$r[s],
                                    lambda = lambdas$M[s])) -> sit
  time_results %<>% rbind(
    data.frame(sparsity = dat$sparsity,
         model = c("IMR", "SoftImpute"),
         time = c(as.numeric(imrt,"sec")[2], as.numeric(sit, "sec")[2])))

    print(time_results)
  }
  saveRDS(time_results, "./notes/saved/sim2_time.rds")
}

time_results %>%
  group_by(model, sparsity) %>%
  summarise_all(c(time_mean=mean,time_sd= sd)) %>%
  arrange(model, sparsity) %>%
  mutate(
    sparsity = 1 - sparsity,
    se   = time_sd / sqrt(30),
    ymin =  pmax(0, time_mean - 1.96 * se) ,
    ymax = (time_mean + 1.96 * se)
  ) -> time_df


ggplot(time_df, aes(x = sparsity, y = time_mean,
                color = model, fill = model, group = model)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_x_continuous("Observation rate", labels = percent_format(accuracy = 1),
                     breaks = seq(0.02, 0.3, 0.04)) +
  scale_y_continuous("Time (s)") +
  scale_color_manual(values = okabe_ito) +
  scale_fill_manual(values  = okabe_ito) +
  labs(linetype = "Model", color = "Model", fill = "Model"
       #title =  "Fit time vs. sparsity (mean ± 95% CI)"
       ) +
  #theme_classic(base_size = 11) +
  theme_bw()+
  theme(legend.position = "top",
        legend.justification = "left",
        strip.text = element_text(face = "bold"))->sim2.g2;sim2.g2

ggsave("./notes/saved/sim2_plot2.png", sim2.g2, width = 320/25.4, height = 150/25.4, dpi = 600)


#
#     sim2.res %>%
#     dplyr::mutate(sparsity = round(sparsity, 2)) %>%
#     dplyr::group_by(model, sparsity, n) %>%
#     dplyr::summarize_all(compute.mean.sd) %>%
#     dplyr::ungroup() %>%
#     arrange(sparsity) %>%
#     mutate(sparsity2 = sparsity) %>%
#     mutate(sparsity = paste0(round(100-sparsity*100),"%")) %>%
#     pivot_longer(c("beta", "M", "theta", "test", "rank", "sparsity"), names_to = c("metric")) ->
#     sim2.tab
#
# # saveRDS(sim2.tab, "./notes/saved/tab_sim2.rds")
# metric_labels <- c(
#   beta = "RMSE($\\beta$)",
#   M    = "RMSE($M$)",
#   theta    = "RMSE($\\Theta$)",
#   test = "Test error",
#   rank       = "Rank($\\Theta$)",
#   sparsity = "Observation rate ($\\%$)"
# )
#
# # method_order <- c(
# #   "IMR", "SoftImpute"
# # )
#
# wide <- sim2.tab %>%
#   mutate(metric = recode(metric, !!!metric_labels)) %>%
#   # mutate(sparsity = round(100-sparsity*100)) %>%
#   # mutate(sparsity = factor(sparsity, levels = seq(30,2, -2),
#   #                          labels=paste0(seq(30,2, -2),"%"))) %>%
#   # mutate(model = factor(model, levels = method_order)) %>%
#   arrange(model, sparsity2) %>%
#   mutate(sparsity3 = factor(sparsity2, levels = seq(.7, .98, .02))) %>%
#   select(model, n, sparsity2, sparsity3, metric, value) %>%
#   pivot_wider(id_cols = c(model, n, sparsity2, sparsity3),
#               names_from = metric, values_from = value) %>%
#   arrange(model, n, sparsity3)
#
# panel_counts <- wide %>% count(model, name = "rows")
# panel_index  <- stats::setNames(panel_counts$rows, panel_counts$model)
#                                 # paste0("$", round(panel_counts$sparsity*100), "\\%$ Missing"))
# require(kableExtra)
# # Build table
# kbl(
#   wide %>% select(-model, -sparsity2, -sparsity3, -n),
#   booktabs  = TRUE,
#   escape    = FALSE,  # keep LaTeX math in headers
#   # col.names = c("Method", setdiff(names(wide), c("n", "method"))),
#   align     = c("l", rep("r", ncol(wide) - 2)),
#   caption   = paste("<span style='color:#000'>",
#                     "RMSE, test error, estimated ranks, and their standard errors (in parenthes",
#                     "under model $\\Theta=X\\beta+M$ and with rank(M)=10, number of covariates = 20,",
#                     "true rank of $\\Theta=30$, 20% observation rate of $\\Theta$, (n,m)=(600,600),(800,800),",
#                     "and 50% of $\\beta$ is 0, at random.</span>")
# ) |>
#   pack_rows(index = panel_index) |>
#   kable_styling(latex_options = c("hold_position", "striped"), font_size = 12) |>
#   row_spec(0, color = "#000") |>
#   row_spec(1:nrow(wide), color = "#000")







#--------------------------
fit <- fit.nlrr$fit
sim1_res(dat, fit2)

fit2 <- IMR::imr.fit(dat$fit_data$Y_full, dat$fit_data$X$Q,
                     lambda_M = fit.nlrr$lambda_M, r = fit.nlrr$rank_M,
                     lambda_beta = fit.nlrr$lambda_beta, warm_start = fit,
                     trace = T)



fit.nlrr$v <- fitsi$fit$v


fit0 <- IMR::imr.fit_no_low_rank(dat$fit_data$train, dat$fit_data$X$Q,
                                 lambda_beta=0.4, intercept_row = T,
                                 intercept_col = T)
fit0$gamma0


partial_crossprod(matrix(1,nrow(dat$Y),1),
                  matrix(fit0$gamma0, nrow=1),
                  dat$fit_data$valid@i, dat$fit_data$valid@p)
quick_camc_simu_res(dat, fit32$fit)

ttrain <- tvalid <- fit.nlrr$fit$resid
ttrain[dat$fit_data$train==0] = 0
tvalid[dat$fit_data$valid==0] = 0
fitc <- simpute.cv(ttrain, tvalid, fit.nlrr$fit$resid, trace=T, tol=15)
fitc$fit$beta <- fit.nlrr$fit$beta
fitc$fit$beta0 <- fit.nlrr$fit$beta0
sim1_res(dat, fitc$fit)


quick_camc_simu_res(dat, fitsi$fit)

fitmm <- IMR::imr.cv_M(ttrain, tvalid, Y_full = fit.nlrr$fit$resid, trace = T)
fitmm$fit$beta <- fit.nlrr$fit$beta
sim1_res(dat, fitmm$fit)

sim1_res(dat, fit32$fit)
sim1_res(dat, o)



quick_camc_simu_res(dat, fitmmci$fit, T)
