% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{get_tie_starts_and_ends}
\alias{get_tie_starts_and_ends}
\title{Identify indices of ties observations in a vector of
ordered outcome values. Not meant to be called directly
by the user.}
\usage{
get_tie_starts_and_ends(ordered_outcomes)
}
\arguments{
\item{ordered_outcomes}{A vector of numeric values arranged
in nondecreasing order.}
}
\value{
A list containing two elements, \code{starts} and
    \code{ends}. For each duplicated value in the input in
    order of appearance, vector \code{starts} gives the
    first index containing the value, and vector \code{ends}
    gives the last (contiguous) index containing the value.
}
\description{
Identify indices of ties observations in a vector of
ordered outcome values. Not meant to be called directly
by the user.
}
