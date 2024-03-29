% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anno.R
\name{spl_mapFromTranscripts}
\alias{spl_mapFromTranscripts}
\title{Map From a transcript to the genome, splitting elements by exons}
\usage{
spl_mapFromTranscripts(trspacegr, exons_grl)
}
\arguments{
\item{trspacegr}{GRanges; an object in transcript space, to be mapped back
to the genome}

\item{exons_grl}{exonsgrl; exons making up the space element is to be mapped
from.}
}
\value{
a granges object containing 1 or more element for each
transcript space range, in genome space, corresponding to pieces
of each element split by exon boundaries
}
\description{
Map From a transcript to the genome, splitting elements by exons
}
\author{
Dermot Harnett, \email{dermot.p.harnett@gmail.com}
}
\keyword{Ribostan}
