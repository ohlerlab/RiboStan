% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kl_div.R
\name{get_kl_df}
\alias{get_kl_df}
\title{Calculate KL-divergence from a data frame of rust expected and actual scores}
\usage{
get_kl_df(metacodondf, anno)
}
\arguments{
\item{metacodondf}{a data frame with metacodon info}

\item{anno}{An annotation object}
}
\value{
A data frame with expected and actual rust score per position in
window, codon, readlength
}
\description{
Calculate KL-divergence from a data frame of rust expected and actual scores
}
\examples{
data(chr22_anno)
data(rpfs)
data(offsets_df)
covgrs <- list(sample1 = rpfs)
# note this doesn't work that well on a small subset
kl_df <- get_kl_df(metacodondf, chr22_anno)
}
\author{
Dermot Harnett, \email{dermot.p.harnett@gmail.com}
}
\keyword{Ribostan}
