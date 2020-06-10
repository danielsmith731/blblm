#' @import purrr
#' @import stats
#' @import readr
#' @import parallel
#' @import tidyverse
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' Linear Regression BLB Original version
#' @param formula regression formula
#' @param data a list of files, or a vector
#' @param m prefered size of data seperation
#' @param B number of bootstraps
#' @export
blblm <- function(formula, data, m = 10, B = 5000) {
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#blblm parallelized for realz homie
#' Title
#'
#' @param formula regression formula
#' @param data a list of files, or a vector
#' @param m prefered size of data seperation
#' @param B number of bootstraps
#' @param cluster parallel cluster
#' @return
#' @export
blblm_par <-function(formula, data, m = 10, B = 5000, cluster){
  data_list <- split_data(data, m)
  estimates <-parLapply(cluster,data_list,
                        function(formula = formula, data = data, n = nrow(data), B=B){
                          lm_each_subsample(formula = formula,data=data,n=n,B=B)
                        })
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

#' split data into m parts of approximated equal sizes
#' @param data a list of files, or a vector
#' @param m prefered size of data seperation
split_data <- function(data, m) {
  if(class(data) == "character"){
    file.path(folder_directory, list.files(folder_directory, pattern = "csv$")) %>%
      map(read_csv)
  }
  else{
    idx <- sample.int(m, nrow(data), replace = TRUE)
    data %>% split(idx)
  }
}


#' LR: compute the estimates
#' @param formula regression formula
#' @param data data, just data
#' @param n number of random vectors
#' @param B number of bootstraps
#' @export
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute the linear regression estimates for a blb dataset
#' @param formula regression formula
#' @param data data, just data
#' @param n number of random vectors
#' @export
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
#' @param formula regression formula
#' @param data data, just data, probably a vector
#' @param freqs weights of model fit
#' @export
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}



blbcoef <- function(fit) {
  coef(fit)
}

glblbcoef <- function(fit) {
  coef(fit)
}

#' compute sigma from fit
#' @param fit stylish clothes you wear for special occasions
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' Generalized Linear Model with BLB
#' @param formula regression formula
#' @param data a list of files, or a vector
#' @param m prefered size of data seperation
#' @param B number of bootstraps
#' @param family the glm family
#' @export
blbglm <- function(formula, data, m = 10, B = 5000, family) {
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ glm_each_subsample(formula = formula, data = ., n = nrow(.), B = B, family))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbglm"
  invisible(res)
}

#' same same but different
#' @param formula regression formula
#' @param data data, just data
#' @param n number of random vectors
#' @param B number of bootstraps
#' @param family family in glm
#' @export
glm_each_subsample <- function(formula, data, n, B, family) {
  replicate(B, glm_each_boot(formula, data, n, family), simplify = FALSE)
}

#' same same but different
#' @param formula regression formula
#' @param data data, just data
#' @param n number of random vectors
#' @param family family in glm
glm_each_boot <- function(formula, data, n, family) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs, family)
}

#' same but different
#' @param formula regression formula
#' @param data data, just data
#' @param freqs weights
#' @param family family in glm
glm1 <- function(formula, data, freqs, family) {
  # drop the original closure of formula,
  # otherwise the formula will pick wrong variables from a parent scope.
  environment(formula) <- environment()
  fit <- glm(formula, data, family = family, weights = freqs)
  list(coef = glblbcoef(fit), sigma = sigma(fit))
}

#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' Regression Model Sigma: BLB
#' @param object Model
#' @param confidence True or False
#' @param level confidence percentage
#' @param ... extra stuff
#' @export sigma.blblm
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' confidence interval
#' @param object Model
#' @param parm TRUE or FALSE
#' @param level confidence percentage
#' @param ... extra stuff
#' @export confint.blblm
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(fit$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}
#' predict
#' @param object Model
#' @param new_data new studff
#' @param confidence TRUE or FALSE
#' @param level confidence percentage
#' @param ... extra conditions
#' @export predict.blblm
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
