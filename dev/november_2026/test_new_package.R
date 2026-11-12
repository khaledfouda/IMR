source("./article/imr-reproducible-code/code/simulation/helper.R")


#devtools::build()
devtools::document()
devtools::load_all(recompoile=TRUE)




simd <- generate_simulated_data(300,200,r=4,p=3,q=3,sparsity=0.8,
                                sparsity_beta=0, sparsity_gamma=0,
                                snr = Inf, mechanism="mcar",
                                identified = TRUE,
                                row_intercept = TRUE,
                                col_intercept = TRUE,
                                seed = 2025)


