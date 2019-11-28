# NOTE: https://stats.stackexchange.com/questions/46532/cox-baseline-hazard
context("Hazard estimation with HAL based on penalized Cox.")
set.seed(74296)
library(glmnet)
library(survival)

# create survival data structures
data(kidney)
y_surv <- Surv(kidney$time, kidney$status)
x_surv <- kidney[, c("age", "sex", "disease", "frail")]
x_surv$disease <- as.numeric(x_surv$disease)
x_surv <- as.matrix(x_surv)

# fit CV-lasso and standard lasso with Cox penalty
cv_fit <- cv.glmnet(x = x_surv, y = y_surv, family = "cox", maxit = 1000)
fit <- glmnet(x = x_surv, y = y_surv, family = "cox", maxit = 1000)

# fit Cox model for baseline hazard (actually gives cumulative baseline hazard)
cph <- coxph(Surv(time, status) ~ age + sex + disease + frail, kidney,
             method = "breslow")
lambda0_cum <- basehaz(cph, centered = FALSE)

# convert from cumulative baseline hazard to standard baseline hazard
lambda0 <- c(lambda0_cum$hazard[1], diff(lambda0_cum$hazard))
lambda0_aug <- rep(lambda0, times = table(sort(kidney$time)))

# try with hal9001 instead of glmnet
hal_cph <- fit_hal(X = x_surv, Y = y_surv, fit_type = "glmnet",
                   family = "cox", max_degree = NULL, return_lasso = TRUE,
                   yolo = FALSE)

# predict from the Cox model
glmnet_pred <- as.numeric(predict(cv_fit, x_surv))
hal_pred <- predict(hal_cph, new_data = x_surv)
