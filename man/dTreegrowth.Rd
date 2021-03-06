% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools2sIPM.R
\name{dTreegrowth}
\alias{dTreegrowth}
\title{dTreegrowth: distribution of dbhs t + 1 as a function
of initial dbh}
\usage{
dTreegrowth(dbhref, dbhtar, model, param, sigma, gridsize, L, U, It = 1,
  Elt = FALSE)
}
\arguments{
\item{dbhref}{tree initial dbh}

\item{dbhtar}{tree target dbh}

\item{model}{fitted model for tree growth and survival}

\item{gridsize}{relavent IPM gridsize}

\item{L}{lower limit of the IPM}

\item{It}{Iteration for numerical sensitivity calculation.
A vector of values, corresponding to parameter dbh/hght and of the same
length, with with to interate the vital rate function. Defaults to 1,
which means no iteration will be done. The function ItMat can be
use to build iterations. see ?ItMat.}

\item{Elt}{only return the parameter relavent for elasticity}
}
\description{
dTreegrowth: distribution of dbhs t + 1 as a function
of initial dbh
}
\author{
Marco D. Visser et al.
}

