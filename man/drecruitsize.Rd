% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools2sIPM.R
\name{drecruitsize}
\alias{drecruitsize}
\title{drecruitsize: density distribution of recruits}
\usage{
drecruitsize(hght, modelparam, It = 1)
}
\arguments{
\item{hght}{seedling height}

\item{modelparam}{a list with the model name}

\item{It}{Iteration for numerical sensitivity calculation. Also see ?ItMat.
A vector of values, corresponding to parameter dbh/hght and of the same
length, with with to interate the vital rate function. Defaults to 1.}
}
\description{
drecruitsize: density distribution of recruits
}
\author{
Marco D. Visser et al.
}

