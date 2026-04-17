source("./article/bixi/data_bixi.R")
source("./article/bixi/helpers_bixi.R")
library(dplyr)
library(magrittr)
dat <- prepare_bixi_data()
X <- dat$X %>% as.matrix()
Z <- dat$Z %>% as.matrix()
Y <- dat$modd$Y
kernels <- generate_similarity_bixi(0.25, "Feb_last", prefix = "1",
                                    train_prefix = "65",
                                    "./article/bixi/data/splits2")
bixi_example <- list(Y=as.matrix(Y), X=X,Z=Z, kernels = kernels, test = as.matrix(dat$test))
usethis::use_data(bixi_example, overwrite = TRUE)
