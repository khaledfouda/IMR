source("./article/Bixi/data_bixi.R")
source("./article/Bixi/helpers_bixi.R")
library(dplyr)
library(magrittr)
dat <- prepare_bixi_data()
X <- dat$X %>% as.matrix()
Z <- dat$Z %>% as.matrix()
Y <- dat$modd$Y
kernels <- generate_similarity_bixi(0.25, "Feb_last", prefix = "1",
                                    train_prefix = "65",
                                    "./article/bixi/data/splits",
                                    return_distance = TRUE)

bixi_example <- list(Y=as.matrix(Y), X=X,Z=Z,
                     test = as.matrix(dat$test),
                     spatial_distance = kernels$distance$spatial,
                     temporal_positions = kernels$distance$temporal)
usethis::use_data(bixi_example, overwrite = TRUE, compress = "xz")

n <- 100; m <- 150
rind <- sample.int(nrow(IMR::bixi_example$Y),n)
cind <- sample.int(ncol(IMR::bixi_example$Y),m)

bixi_sample <- list( Y = IMR::bixi_example$Y[rind,cind],
                     test = IMR::bixi_example$test[rind,cind],
                      X = as.matrix(IMR::bixi_example$X[rind,]),
                      Z = as.matrix(IMR::bixi_example$Z[cind,]),
                      spatial_distance = IMR::bixi_example$spatial_distance[cind,cind],
                      temporal_distance = IMR::bixi_example$temporal_positions[rind,rind])


usethis::use_data(bixi_sample, overwrite = TRUE, compress = "xz")


saveRDS(IMR::bixi_example, "article/bixi_example.rds")
