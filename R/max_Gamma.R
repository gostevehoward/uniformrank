#' Report the largest sensitivity parameter Gamma for which the
#' sensitivity analysis for a uniform general signed rank test rejects
#' the null hypothesis at level alpha., i.e. the largest value of
#' Gamma for which the upper uniform boundary intersects the random
#' walk induced by the observed data.
#'
#' @param outcomes A vector of matched pair differences.
#' @param score_fn A string indicating which score function should be
#'      used.  Must be one of 'sign' (sign score function), 'wsrt'
#'      (Wilcoxon signed rank score function), 'normal_scores'
#'      (normal score function) and 'redescending'.
#' @param alpha Rate of type I error control for the hypothesis test.
#' @param x0 Tuning parameter for the uniform testing boundary, chosen
#'     in the range $[0,1]$.
#' @return The maximum sensitivity parameter value (a single numeric
#'     value no smaller than 1).
#' @details Ties are handled automatically.  As described in greater
#'     detail in the manuscript referenced below, all tied observations
#'     are grouped into a single step of the random walk.
#' @examples
#' data(mercury2)
#' fixed_max_Gamma(mercury2, 'sign')
#' uniform_max_Gamma(mercury2, 'sign')
#' @references Howard, S.R., and Pimentel, S.D. (2020), "The uniform general
#' signed rank test and its design sensitivity." arXiv:1904.08895.
#' @name uniform_max_Gamma
uniform_max_Gamma <- function(outcomes, score_fn = 'sign', alpha = 0.05, x0=1/3) {
  score.out <- get_scores_and_T_by_score(outcomes, score_fns(score_fn))
  scores <- score.out$scores
  Tvec <- score.out$T
  n0 <- floor(x0 * length(Tvec))
  c <- rev(scores)
  sum_csq_n0 <- sum(c[1:n0]^2)

  root_equation <- function(rho) {
    lambda <- sqrt(2 * log(1 / alpha) / (sum_csq_n0 * rho * (1 - rho)))
    boundary <- (
      1 / lambda
      * (log(1 / alpha) + cumsum(log(1 + rho * (exp(c * lambda) - 1)))))
    min(boundary - Tvec)
  }

  compute_Gamma_from_rho(root_equation)
}


#' Report the largest sensitivity parameter Gamma for which the
#' sensitivity analysis for a uniform general signed rank test rejects
#' the null hypothesis at level alpha., i.e. the largest value of
#' Gamma for which the upper uniform boundary intersects the random
#' walk induced by the observed data.
#'
#' @param outcomes A vector of matched pair differences.
#' @param score_fn A string indicating which score function should be
#'      used.  Must be one of 'sign' (sign score function), 'wsrt'
#'      (Wilcoxon signed rank score function), 'normal_scores'
#'      (normal score function) and 'redescending'.
#' @param alpha Rate of type I error control for the hypothesis test.
#' @return The maximum sensitivity parameter value (a single numeric
#'     value no smaller than 1).
#' @examples
#' data(mercury2)
#' fixed_max_Gamma(mercury2, 'sign')
#' uniform_max_Gamma(mercury2, 'sign')
#' @references Howard, S.R., and Pimentel, S.D. (2020), "The uniform general
#' signed rank test and its design sensitivity." arXiv:1904.08895.
#' @name fixed_max_Gamma
fixed_max_Gamma <- function(outcomes, score_fn = 'sign', alpha = 0.05) {
  score.out <- get_scores_and_T_by_score(outcomes, score_fns(score_fn))
  scores <- score.out$scores
  final_T <- tail(score.out$T, 1)
  unscaled_mean <- sum(scores)
  unscaled_var <- sum(scores^2)
  z_factor <- qnorm(1 - alpha)

  root_equation <- function(rho) {
    (rho * unscaled_mean
     + z_factor * sqrt(rho * (1 - rho) * unscaled_var)
     - final_T)
  }
  compute_Gamma_from_rho(root_equation)
}

