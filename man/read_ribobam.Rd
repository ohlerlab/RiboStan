% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ritpms.R
\name{read_ribobam}
\alias{read_ribobam}
\title{Read a bam file containing Ribosomal Footprints}
\usage{
read_ribobam(ribobam, which, strip_seqnames = TRUE)
}
\arguments{
\item{ribobam}{String; full path to html report file.}

\item{which}{Select which reads to include (e.g. those that overlap
with start and stop positions)}

\item{strip_seqnames}{whether the function should remove all text after the
first '|', useful if aligning to gencode fastas#
Defaults to \code{TRUE}}
}
\value{
This function returns a granges object with the read names in the
names slot and a metadata column
denoting readlength.
}
\description{
This function
}
\details{
This function reads in a bam file containing ribosomal footprints.
It uses only reads without splice sites and reads which align to the positive
strand, as it's designed to work on transcriptomic alignments.
}
\seealso{
\code{\link{get_cds_reads}}, \code{\link{get_readlens}}
}
\author{
Dermot Harnett, \email{dermot.p.harnett@gmail.com}
}
\keyword{Ribostan}
