% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsaHits.R
\name{rsaHits}
\alias{rsaHits}
\title{Select hits by RSA}
\usage{
rsaHits(
  dta,
  treatment,
  control,
  normMethod = "PLATE",
  LB,
  UB,
  revHits = FALSE,
  Bonferroni = FALSE,
  outputFile = "RSAhits.csv",
  scoreFile = "RSA_score.csv"
)
}
\arguments{
\item{dta}{synthetic lethal RNAi screen data}

\item{treatment}{the treatment condition in EXPERIMENT_MODIFICATION}

\item{control}{the control condition in EXPERIMENT_MODIFICATION}

\item{normMethod}{normalization methods. If "PLATE", then values are normalized by plate median, otherwise use the
provided control siRNA}

\item{LB}{Low bound}

\item{UB}{up bound}

\item{revHits}{reverse hit picking, default the lower the score the better}

\item{Bonferroni}{conceptually useful when there are different number of siRNAs per gene, default FALSE}

\item{outputFile}{output file name}

\item{scoreFile}{name of the score file to be written under the current folder}
}
\value{
A result file written to the current folder.
\itemize{
  \item Gene_ID,Well_ID,Score: columns from input spreadsheet
  \item LogP: OPI p-value in log10, i.e., -2 means 0.01
  \item OPI_Hit: whether the well is a hit, 1 means yes, 0 means no
  \item #hitWell: number of hit wells for the gene
  \item #totalWell: total number of wells for the gene. If gene A has three wells w1, w2 and w3, and w1 and w2 are hits, #totalWell should be 3, #hitWell should be 2, w1 and w2 should have OPI_Hit set as 1 and w3 should have OPI_Hit set as 0.
  \item OPI_Rank: ranking column to sort all wells for hit picking
  \item Cutoff_Rank: ranking column to sort all wells based on Score in the simple activity-based method
}
Note: a rank value of 999999 means the well is not a hit
}
\description{
Selected hits by redundant siRNA activity method. Here is a wrapper function of RSA 1.8 by Yingyao Zhou.
}
\examples{
data(example_dt)
rsaHits(example_dt, treatment = "treatment", control = "control",
        normMethod = "PLATE", LB = 0.2, UB = 0.8, revHits = FALSE,
        Bonferroni = FALSE, outputFile = "RSAhits.csv")
}
\references{
Koenig, R. et al. A probability-based approach for the analysis of large-scale RNAi screens. Nat Methods 4, 847-849 (2007).
}
