% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zFactor.R
\name{zFactor}
\alias{zFactor}
\title{Calcualte the Z and Z' factor}
\usage{
zFactor(dta, negativeCon, positiveCon, useMean = TRUE)
}
\arguments{
\item{dta}{synthetic lethal RNAi screen data.}

\item{negativeCon}{the negative control used in the WELL_CONTENT_NAME.}

\item{positiveCon}{the positive control used in the WELL_CONTENT_NAME.}

\item{useMean}{use mean to calcualate z factor and z' factor by default; otherwise use median.}
}
\value{
A data.frame contains z factor and z' factor
}
\description{
calcualte the Z and Z' factor for each plate.
}
\examples{
data(example_dt)
res <- zFactor(example_dt, negativeCon = "scrambled control si1", positiveCon = "PLK1 si1")
}
\references{
Zhang J.H., Chung T.D. & Oldenburg K.R. A simple statistical parameter for use in evaluation and validation of high throughput screening assays. J. Biomol. Screen. B, 4 67-73 (1999).
Birmingham,A. et al. (2009) Statistical methods for analysis of high-throughput RNA interference screens. Nat Methods, 6, 569-575.
}
