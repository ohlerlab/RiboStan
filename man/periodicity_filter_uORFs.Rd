% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/periodicity_tests.R
\name{periodicity_filter_uORFs}
\alias{periodicity_filter_uORFs}
\title{Run Multitaper tests on a set of ORFs}
\usage{
periodicity_filter_uORFs(psites, anno, remove = TRUE, n_cores = 1)
}
\arguments{
\item{psites}{GRanges object with psite information}

\item{anno}{annotation object}

\item{remove}{whether to remove non-peridic uORFs}

\item{n_cores}{number of cores to use}
}
\value{
a numeric vector with the spectral coefficient at 0.333... and the pvalue for the test
}
\description{
Run Multitaper tests on a set of ORFs
}
\details{
This function applies a multitaper test to
}
\examples{
data(chr22_anno)
data(rpfs)
data(offsets_df)
psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
filteredanno <- periodicity_filter_uORFs(psites, chr22_anno)
}
\author{
Dermot Harnett, \email{dermot.p.harnett@gmail.com}
}
\keyword{Ribostan}
