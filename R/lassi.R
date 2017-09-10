#' generate sequence of lambdas
#' @export
lambda_seq <- function(lambda_max, lambda_min_ratio = 0.01, nlambda = 100){
  log_seq <- seq(from=0,to=log10(lambda_min_ratio), length=nlambda)
  result <- lambda_max * 10^log_seq
  return(result)
}
