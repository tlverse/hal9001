context("Fits and prediction of classic Super Learner with HAL.")

# simulate data
p <- 5
n <- 1000
x <- as.data.frame(replicate(p, rnorm(n)))
y <- rnorm(n, mean = mean(rowMeans(x)), sd = sd(sin(x[, p])))

test_x <- as.data.frame(replicate(p, rnorm(n)))
test_y <- rnorm(n, mean = mean(rowMeans(x)), sd = sd(sin(x[, p - 1])))

# run HAL by itself
hal <- fit_hal(X = x, Y = y, yolo = FALSE)
pred_hal_train <- predict(hal, newX = x)  # returns a numeric vector
pred_hal_test <- predict(hal, newX = test_x)  # returns a numeric vector

# run SL-classic with glmnet and get predictions
glmnet_sl <- SuperLearner(Y = y, X = x, SL.lib = "SL.glmnet")
pred_glmnet_train <- predict(glmnet_sl, newX = x)  # returns a list
pred_glmnet_test <- predict(glmnet_sl, newX = test_x)  # returns a list

# run SL-classic with glmnet and get predictions
hal_sl <- SuperLearner(Y = y, X = x, SL.lib = "SL.hal9001")
pred_hal_sl_train <- predict(hal_sl, newX = x)  # returns a list
pred_hal_sl_test <- predict(hal_sl, newX = test_x)  # returns a list

# test via equivalence of outputs: HAL vs. SL-HAL
expect_equal(pred_hal_train, as.numeric(pred_hal_sl_train$library.predict))
expect_equal(pred_hal_test, as.numeric(pred_hal_sl_test$library.predict))

# test via equivalence of outputs: HAL vs. SL-glmnet
expect_equal(pred_hal_train, as.numeric(pred_glmnet_train$library.predict))
expect_equal(pred_hal_test, as.numeric(pred_glmnet_test$library.predict))

# test via equivalence of outputs: SL-HAL vs. SL-glmnet
expect_equal(as.numeric(pred_hal_sl_train$library.predict),
             as.numeric(pred_glmnet_train$library.predict)
            )
expect_equal(as.numeric(pred_hal_sl_test$library.predict),
             as.numeric(pred_glmnet_test$library.predict)
            )

