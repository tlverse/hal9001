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

# fit Cox model for baseline hazard (actually gives cumulative baseline hazard)
cph <- coxph(Surv(time, status) ~ age + sex + disease + frail, kidney,
  method = "breslow"
)
lambda0_cum <- basehaz(cph, centered = FALSE)

# convert from cumulative baseline hazard to standard baseline hazard
haz <- exp(diff(lambda0_cum[, 1]) * diff(lambda0_cum[, 2]))
lambda0 <- rep(c(lambda0_cum$hazard[1], diff(lambda0_cum$hazard)),
  times = table(sort(kidney$time))
)

# fit CV-lasso with Cox penalty and predict
cv_coxnet <- cv.glmnet(x = x_surv, y = y_surv, family = "cox")
coxnet_pred <- as.numeric(predict(cv_coxnet, x_surv, type = "response"))

# try with hal9001 instead of glmnet
cv_halcox <- fit_hal(
  X = x_surv, Y = y_surv, fit_type = "glmnet",
  family = "cox", yolo = FALSE
)
halcox_pred <- predict(cv_halcox, new_data = x_surv)



# fit lasso with Cox penalty over a grid of lambda and predict
nocv_coxnet <- glmnet(x = x_surv, y = y_surv, family = "cox", nlambda = 200)
nocv_coxnet_pred <- as.matrix(predict(nocv_coxnet, x_surv, type = "response"))

# fit HAL with Cox penalty over a grid of lambda and predict
nocv_halcox <- fit_hal(
  X = x_surv, Y = y_surv, fit_type = "glmnet", # ,
  family = "cox", cv_select = FALSE, nlambda = 200,
  yolo = FALSE
)
nocv_halcox_pred <- predict(nocv_halcox, new_data = x_surv)
