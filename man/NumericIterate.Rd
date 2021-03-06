% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools2sIPM.R
\name{NumericIterate}
\alias{NumericIterate}
\title{Numerically Iterate size-specific vital rates and retrive
the dominant eigenvalue.}
\usage{
NumericIterate(seedlingclasses, treeclasses, seedlingclasseswidth,
  treeclasseswidth, sizerange, hghtdbhthreshold, dbhthreshold, f.sp, f.sts,
  f.repro, f.sdlhght, f.s2t, f.sg, f.ss, f.tg, f.ts,
  vitalrate = "tree growth", delta = 1.01)
}
\arguments{
\item{vitalrate}{which vital rate to iterate? See ?ItMat.}

\item{delta}{the factor of iteration
All parameters are the same as in the constructIPM function.}
}
\description{
Numerically Iterate size-specific vital rates and retrive
the dominant eigenvalue.
}

