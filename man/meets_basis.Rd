% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{meets_basis}
\alias{meets_basis}
\title{Compute Values of Basis Functions}
\usage{
meets_basis(X, row_num, cols, cutoffs, orders)
}
\arguments{
\item{X}{The design matrix, containing the original data.}

\item{row_num}{Numeri for  a row index over which to evaluate.}

\item{cols}{Numeric for the column indices of the basis function.}

\item{cutoffs}{Numeric providing thresholds.}

\item{orders}{Numeric providing smoothness orders}
}
\description{
Computes and returns the indicator value for the basis described by
cols and cutoffs for a given row of X
}
