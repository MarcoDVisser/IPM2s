% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools2sIPM.R
\name{Rkernel}
\alias{Rkernel}
\title{Rkernel: Recruits kernel}
\usage{
Rkernel(dbhref, hghttar, models, gridsize, treegridsize, smodels, tImat)
}
\arguments{
\item{dbhref}{mother tree diameter}

\item{hghttar}{recruit hght}

\item{models}{list of models or parameters for seed 
to seedling transition, reproduction, seed production, 
and inititail recruit height}

\item{gridsize}{gridsize of the seedling IPM}

\item{smodels}{Tree survival models}

\item{tImat}{Tree iteration matrix, which in this case
is expected to be a vector with values for seed production, reproduction,
recruitsize  tree growth &  tree survival}
}
\description{
Function returns recruitment (N individuals) between 
dbhref and hghttar
}
\author{
Marco D. Visser et al.
}

