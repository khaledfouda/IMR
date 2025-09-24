
fit_MA25_movielens <- function(input_tag = "_c_10_98_",
                               seed = 2025){

  if(is.numeric(seed)) set.seed(seed)

  load("./notes/movielens/data/Movie_X.Rdata") #X
  X <- X[,1:4]
  if(input_tag == ""){
    load("notes/movielens/data/Movie_Y.Rdata",verbose = T)
  }else
    Y <-     readRDS(paste0("./notes/movielens/data/Movie_Y",input_tag,".Rdata"))
  # query <- readRDS(paste0("./notes/movielens/data/ml-1m/Movie_Q",input_tag,".Rdata"))

  # load original MA functions
  source('./other_models/Ma25.R')
  #============
  file_location <- paste0("./notes/movielens/data/saved_models/Ma_fit",input_tag,".rds")
  M <- MA25.fit(Y, X, save_to_file = T, file_location = file_location)
  return(M)

    #
  # # Initialize the solver object (MAR with logistic missing model)
  # M <- MC_alt_LS_rank_solver(Y=Y, X=X, intercept_val=1,
  #                            max_iter = 30L, .tol = 1e-6, missing_model = 'logistic')
  # rm(X, Y)
  #
  # Gender <- 'G'
  # Age <- paste0('A', 1:3)
  # M$.make_obs()
  #
  # #===========================================================
  # # Table 6: The fitting result for the logistic model for π_i
  # #===========================================================
  #
  # # Estimate the missing model
  # M$estimate_pi()
  # R <- summary(M$pi_est$fit)
  # Tab6 <- R$coefficients
  # rownames(Tab6) <- c('(Intercept)', Gender, Age, paste(Gender, Age, sep='*'))
  # print(Tab6)
  # #===========================================================
  # # Table 7: RMSE of different models in training and test set
  # #===========================================================
  # # Include only main effects for model fitting (best performing based on the article)
  # M$set_fitting_cov(index = c(1:5))
  # M$summary()
  # # Rank estimation
  # M$rank_estimation(C_h=0.2, delta_h=0.1, penalty_names = c('h'), r_bar=8)
  # rhat <- M$rank_est$est['h'] # get the rank estimation result
  # # Fitting with the estimated rank
  # M$fitting(with_rs = rhat)
  # saveRDS(M, paste0("data/saved_models_1m/Ma_fit",input_tag,".rds"))
}
