#' Gencode v32 annotation for human chromosome 22
#'
#' The gencode v32 annotation, downloaded from Annotationhub()
#' and then processed with load_annotation, including uORFs
#'
#' @format An annotation object as created by load_annotation:
#' @source \url{AnnotationHub gencode v32 annotation}
"chr22_anno"

#' A data frame containing metacodon data
#' @format A tibble: with 23,058 rows and 8 columns
#' \describe{
#'   \item{sample}{the sample used to calculate the metacodons}
#'   \item{codon}{the codon used to calculate the metacodon profile}
#'   \item{position}{the position relative tot he codon}
#'   \item{ro_cl}{The RUST score at that position}
#'   \item{re_c}{The expected RUST score at that position}
#'   \item{count}{the mean read 5' end count at that position}
#' }
#'
#' @source \url{AnnotationHub gencode v32 annotation}
"metacodondf"

#' GRanges object with 63346 ranges and 1 metadata column:
#'
#' @format GRanges object with 63346 ranges and 1 metadata column:
#' \describe{
#'   \item{readlen}{the readlength of the footprint}
#' }
#' @source \url{Simulated data}
"rpfs"

#' A Dataframe containing CDS-occupancy derived offsets
#' for the RPF data in rpfs
#'
#' @format a  tibble with 12 rows and 3 columns
#' \describe{
#'   \item{readlen}{the readlength to which the offset shoudl be applied}
#'   \item{phase}{the phase to which the offset shoudl be applied}
#'   \item{offset}{the offset which should be applied}
#' }
#' @source \url{AnnotationHub gencode v32 annotation}
"offsets_df"
