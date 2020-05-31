#' Identify indices of ties observations in a vector of
#' ordered outcome values. Not meant to be called directly
#' by the user.
#'
#' @param ordered_outcomes A vector of numeric values arranged
#'     in nondecreasing order.
#' @return A list containing two elements, \code{starts} and
#'     \code{ends}. For each duplicated value in the input in
#'     order of appearance, vector \code{starts} gives the
#'     first index containing the value, and vector \code{ends}
#'     gives the last (contiguous) index containing the value.
#' @name get_tie_starts_and_ends
get_tie_starts_and_ends <- function(ordered_outcomes) {
    ties_pattern <- duplicated(abs(ordered_outcomes))
    search_strings <- paste(as.numeric(ties_pattern),
                            c(ties_pattern[-1], 0), sep = '')
    starts <- which(search_strings == '01')
    ends <- which(search_strings == '10')
    stopifnot(length(starts) == length(ends))

    list(starts=starts, ends=ends)
}


#' Generate the score for each pair difference and the observed
#' value of the random walk T using the raw pair differences and
#' the score function. Not meant to be called directly by users.
#'
#' @param outcomes A vector of matched pair differences.
#' @param score_fn A string indicating which score function should be
#'      used.  Must be one of 'sign' (sign score function), 'wsrt'
#'      (Wilcoxon signed rank score function), 'normal_scores'
#'      (normal score function) and 'redescending'.
#' @return A list containing two elements, \code{scores} and
#'     \code{T}. Vector \code{scores} gives the score for each
#'     pair difference in non-decreasing order. Vector \code{T}
#'     gives the corresponding random walk for the observed data.
#' @name get_scores_and_T_by_score
get_scores_and_T_by_score <- function(outcomes, score_fn) {
  ordered_outcomes <- outcomes[order(abs(outcomes))]
  signs <- (ordered_outcomes > 0)
  scores <- get_scores(score_fn, length(outcomes))
  ties <- get_tie_starts_and_ends(ordered_outcomes)

  replace_scores_with_means <- function(start_index, end_index) {
    scores[start_index:end_index] <<- mean(scores[start_index:end_index])
  }
  mapply(replace_scores_with_means, ties$starts, ties$ends)

  score.order <- order(scores)
  scores <- scores[score.order]
  signs <- signs[score.order]

  Tvec <- cumsum(rev(scores * signs))

  flatten_T_within_tie_group <- function(start_index, end_index) {
    T_start_index <- length(outcomes) - end_index + 1
    T_end_index <- length(outcomes) - start_index + 1
    Tvec[T_start_index:(T_end_index - 1)] <<-
      if (T_start_index == 1) 0 else Tvec[T_start_index - 1]
  }
  mapply(flatten_T_within_tie_group, match(ties$starts, score.order),
         match(ties$ends, score.order))

  list(scores=scores, T=Tvec)
}


#' Given the name of a score function as input, return the score
#' function itself. Not meant to be called directly by users.
#' @param scorename A string indicating which score function should be
#'      used.  Must be one of 'sign' (sign score function), 'wsrt'
#'      (Wilcoxon signed rank score function), 'normal_scores'
#'      (normal score function) and 'redescending'.
#' @return A score function taking a numeric vector of elements in $[0,1]$
#'     and returning another nonnegative numeric vector of the same length.
#' @name score_fns
score_fns <- function(scorename){
    if(scorename == 'sign') return(function(q) rep(1, length(q)))
    if(scorename == 'wsrt') return(function(q) q)
    if(scorename == 'normal_scores') return(function(q) ifelse((1 + q) / 2 >= 1, 0, qnorm((1 + q) / 2)))
    if(scorename == 'redescending') return(get_u_statistic_score_fn(20, 12, 19))
    stop('scorename must be one of sign, wsrt, normal_scores, and redescending')
}

#' The approximate score associated with the u-statistic of Rosenbaum (2011),
#' which counts, for all subsets of size \code{m} of the set of matched pairs,
#' the number of positive differences in the pairs with the \code{m_lower} to
#' \code{m_upper} smallest absolute differences, inclusive. Not meant to be
#' called directly by users.
#' @param m Size of subsets to consider.
#' @param m_lower Lower boundary in the within-subset rank-ordering for pairs to
#'     consider.  Must be no larger than \code{m}.
#' @param m_upper Upper boundary in the within-subset rank-ordering for pairs to
#'     consider.  Must be no larger than \code{m} and no smaller than
#'     $\code{m_lower}.
#' @return A function mapping numeric values in $[0,1]$ to nonnegative numeric
#'     values.
#' @references Rosenbaum, P. R. (2011). A new u-statistic with superior design
#'     sensitivity in matched observational studies. \emph{Biometrics}, 67(3),
#'     1017-1027.
#' @name get_u_statistic_score_fn
get_u_statistic_score_fn <- function(m, m_lower, m_upper) {
    l <- m_lower:m_upper
    Vectorize(function(q) {
        sum(l * choose(m, l) * q^(l-1) * (1 - q)^(m - l) / m)
    })
}

#' Given a score function and a sample size n, returns scores for points
#' spaced equally at 1/(n+1) intervals within the unit interval.  Not meant
#' to be called directly by users.
#' @param score_fn A score function taking a numeric vector of elements in
#'     $[0,1]$ and returning another nonnegative numeric vector of the same
#'     length.
#' @param num_units The number of matched pairs in the study, an integer.
#' @return A vector of length \code{num_units} given the scores for
#'     equispaced points within the unit interval.
#' @name get_scores
get_scores <- function(score_fn, num_units) {
  score_fn((1:num_units) / (num_units + 1))
}

#' Given the maximal probability of treatment within a matched pair in
#' a sensitivity analysis, compute the associated sensitivity parameter
#' Gamma.  Not meant to be called directly by users.
#' @param root_equation A function of treatment probability rho that
#'     achieves value zero for a situation of interest.
#' @return The value of Gamma (a single numeric value) associated with
#'     the value of rho that solves the input root equation.
#' @name compute_Gamma_from_rho
compute_Gamma_from_rho <- function(root_equation) {
  if (root_equation(1/2) > 0) {
    0
  } else {
    rho <- uniroot(root_equation, c(1/2, 1 - 1e-5))$root
    rho / (1 - rho)
  }
}


