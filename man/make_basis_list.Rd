% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{make_basis_list}
\alias{make_basis_list}
\title{Sort Basis Functions}
\usage{
make_basis_list(X_sub, cols, order_map)
}
\arguments{
\item{X_sub}{A subset of the columns of X, the original design matrix.}

\item{cols}{An index of the columns that were reduced to by sub-setting.}

\item{order_map}{A vector with length the original unsubsetted matrix X which specifies the smoothness of the function in each covariate.}
}
\description{
Build a sorted list of unique basis functions based on columns, where each
basis function is a list
}
\details{
Note that sorting of columns is performed such that the basis order
equals cols.length() and each basis function is a list(cols, cutoffs).
}
