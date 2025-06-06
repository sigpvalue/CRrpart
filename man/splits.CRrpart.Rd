% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CRrpart.R
\name{splits.CRrpart}
\alias{splits.CRrpart}
\title{Find the best split  (function not to be called)}
\usage{
splits.CRrpart(x., y., x.var.def, minbuck, failcode, cencode, iter)
}
\arguments{
\item{x.}{Predictor variables}

\item{y.}{Response variables}

\item{x.var.def}{identify categorical predictor variable}

\item{minbuck}{Minimum number of observations needed in each daughter node}

\item{failcode}{Code of fstatus corresponding to the event of interest}

\item{cencode}{Code of fstatus corresponding to censored observations}

\item{iter}{Number of iterations for p-value calculations. (default = 2000)}
}
\description{
Find the best split among potential splits
}
