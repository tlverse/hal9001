context("Unit test for procedures relating to undersmooth HAL.")
# library(microbenchmark)
library(hal9001)
library(glmnet)
library(here)
source(paste0(here(), "/R/hal_undersmooth.R"))

# library(devtools)
# load_all(here())

set.seed(385971)

# simulate data
n <- 100
p <- 3
x <- matrix(rnorm(n * p), n, p)
y <- x[, 1] * sin(x[, 2]) + rnorm(n, mean = 0, sd = 0.2)

# fit undersmoothed HAL
uhal_fit <- fit_uhal(X=x, Y=y, family="gaussian")
uhal_fit$lambda






# # initialize the undersmoothing procedure
# hal_init <- undersmooth_init(X=x, Y=y,
#                              Nlam = 20,
#                              formula = NULL,
#                              X_unpenalized = NULL,
#                              max_degree = ifelse(ncol(x) >= 20, 2, 3),
#                              smoothness_orders = 0,
#                              num_knots = num_knots_generator(
#                                max_degree = ifelse(ncol(x) >= 20, 2, 3),
#                                smoothness_orders = 0,
#                                base_num_knots_0 = 200,
#                                base_num_knots_1 = 50
#                              ),
#                              reduce_basis = 1 / sqrt(length(y)),
#                              family = c("gaussian", "binomial", "poisson", "cox"),
#                              lambda = NULL,
#                              id = NULL,
#                              offset = NULL,
#                              fit_control = list(
#                                cv_select = TRUE,
#                                n_folds = 10,
#                                foldid = NULL,
#                                use_min = TRUE,
#                                lambda.min.ratio = 1e-4,
#                                prediction_bounds = "default"
#                              ),
#                              basis_list = NULL,
#                              return_lasso = TRUE,
#                              return_x_basis = TRUE,
#                              yolo = FALSE)
#
# # do undersmoothed HAL
# hal_under <- undersmooth_hal(X=x,
#                              Y=y,
#                              fit_init=hal_init$fit_init,
#                              basis_mat=hal_init$basis_mat,
#                              Nlam = 20,
#                              family = "gaussian")
#
# hal_under$lambda_init
# hal_under$lambda_under
# hal_under$spec_under
