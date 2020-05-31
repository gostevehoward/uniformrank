#' Report the p-value from a uniform general signed rank test, i.e. the
#' smallest value of Type I error control rate alpha for which the
#' upper uniform boundary intersects the random walk induced by the
#' observed data.
#'
#' @param outcomes A vector of matched pair differences.
#' @param score_fn A string indicating which score function should be
#'      used.  Must be one of 'sign' (sign score function), 'wsrt'
#'      (Wilcoxon signed rank score function), 'normal_scores'
#'      (normal score function) and 'redescending'.
#' @param Gamma Sensitivity parameter (>= 1) describing the degree
#'     of unmeasured confounding present.  The default case Gamma = 1
#'     corresponds to no unmeasured confounding.
#' @param x0 Tuning parameter for the uniform testing boundary, chosen
#'     in the range $[0,1]$.
#' @param tol A minimum tolerance for distinguishing p-values from 0
#'     and 1.  P-values found to be below \code{tol} or above 1 -
#'     \code{tol} will be reported as 0 and 1 respectively.
#' @return A p-value.
#' @details Ties are handled automatically.  As described in greater
#'     detail in the manuscript referenced below, all tied observations
#'     are grouped into a single step of the random walk.
#' @examples
#' data(mercury2)
#' uniform_pval(mercury2, 'sign') #test assuming no unobserved confounding
#' uniform_pval(mercury2, 'sign', Gamma = 5) #with large amount of confounding
#' fixed_pval(mercury2, 'sign', Gamma = 5) #fixed test is more sensitive
#' @references Howard, S.R., and Pimentel, S.D. (2020), "The uniform general
#' signed rank test and its design sensitivity." arXiv:1904.08895.
#' @name uniform_pval
uniform_pval <- function(outcomes, score_fn = 'sign', Gamma = 1, x0 = 1/3,
                         tol = 1e-15){
  score.out <- get_scores_and_T_by_score(outcomes, score_fns(score_fn))
  scores <- score.out$scores
  Tvec <- score.out$T
  n0 <- floor(x0 * length(Tvec))
  c <- rev(scores)
  sum_csq_n0 <- sum(c[1:n0]^2)
  rho <- Gamma/(1 + Gamma)

  root_equation <- function(alpha) {
    lambda <- sqrt(2 * log(1 / alpha) / (sum_csq_n0 * rho * (1 - rho)))
    boundary <- (
      1 / lambda
      * (log(1 / alpha) + cumsum(log(1 + rho * (exp(c * lambda) - 1)))))
    min(boundary - Tvec)
  }
  if(root_equation(tol) < 0 && root_equation(1-tol) < 0) return(0)
  if(root_equation(tol) > 0 && root_equation(1-tol) > 0) return(1 )
  return(uniroot(root_equation, c(tol, 1 - tol))$root)
}


#' Report the p-value from a general signed rank test, or from an
#'     associated sensitivity analysis.
#'
#' @param outcomes A vector of matched pair differences.
#' @param score_fn A string indicating which score function should be
#'      used.  Must be one of 'sign' (sign score function), 'wsrt'
#'      (Wilcoxon signed rank score function), 'normal_scores'
#'      (normal score function) and 'redescending'.
#' @param Gamma Sensitivity parameter (>=1) describing the degree
#'     of unmeasured confounding present.  The default case Gamma = 1
#'     corresponds to no unmeasured confounding.
#' @param tol A minimum tolerance for distinguishing p-values from 0
#'     and 1.  P-values found to be below \code{tol} or above 1 -
#'     \code{tol} will be reported as 0 and 1 respectively.
#' @return A p-value.
#' @examples
#' data(mercury2)
#' uniform_pval(mercury2, 'sign') #test assuming no unobserved confounding
#' uniform_pval(mercury2, 'sign', Gamma = 5) #with large amount of confounding
#' fixed_pval(mercury2, 'sign', Gamma = 5) #fixed test is more sensitive
#' @references Howard, S.R., and Pimentel, S.D. (2020), "The uniform general
#' signed rank test and its design sensitivity." arXiv:1904.08895.
#' @name fixed_pval
fixed_pval <- function(outcomes, score_fn = 'sign', Gamma = 1, tol = 1e-15){
  score.out <- get_scores_and_T_by_score(outcomes, score_fns(score_fn))
  final_T <- tail(score.out$T, 1)
  unscaled_mean <- sum(score.out$scores)
  unscaled_var <- sum(score.out$scores^2)
  rho <- Gamma/(1 + Gamma)

  root_equation <- function(alpha) {
    (rho * unscaled_mean
     + qnorm(1 - alpha) * sqrt(rho * (1 - rho) * unscaled_var)
     - final_T)
  }
  if(root_equation(tol) < 0 && root_equation(1-tol) < 0) return(0)
  if(root_equation(tol) > 0 && root_equation(1-tol) > 0) return(1 )
  return(uniroot(root_equation, c(tol, 1 - tol))$root)
}
