#' Add 'correct offset' column to a set of reads
#'
#' Given rpfs, and named vectors of starts and stops, this adds the 'correct'
#' offset to the reads to use with determining p-site offsets
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpfs GRanges object with RPFs
#' @param anno annotation object
#'
#' @return the rpf granges but wth a 'cor_offset' column


add_cor_offset <- function(rpfs, anno) {
  cdsstarts <- anno$cdsstarts%>%unlist
  cds_prestop_st <- anno$cds_prestop_st%>%unlist

  rpfs <- rpfs %>%
    as("GRanges") %>%
    addphase(cdsstarts) %>%
    {
      .$startoffset <- cdsstarts[as.vector(seqnames(.))] + .$phase - start(.)
      .$endoffset <- cds_prestop_st[as.vector(seqnames(.))] + .$phase - start(.)
      .$cor_offset <- ifelse(
        (.$startoffset < (.$readlen - 5)) & (.$startoffset > 5),
        .$startoffset,
        ifelse(
          (.$endoffset < (.$readlen - 5)) & (.$endoffset > 5),
          .$endoffset,
          NA
        )
      )
      .
    }
  rpfs
}



#' Get Offsets by maximizing p-sites within the CDS
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpfs GRanges object with RPFs
#' @param anno annotation object
#'
#' @return a data frame with the best offset per read/phase

# TODO document the columns output by this function
# TODO maybe also shortedn this by removing arguments?
# TODO remove the hardcoded numbers

get_incl_max_offsets <- function(rpfs, anno)

{
  stopifnot("readlen" %in% colnames(mcols(rpfs)))
  #
  rpfs <- rpfs %>% add_cor_offset(anno)
  #
  rpfs$readlen <- as.numeric(rpfs$readlen)
  if (!"score" %in% colnames(mcols(rpfs))) rpfs$score <- 1
  # data frame with numbers of reads gained at a specific offset
  gaindf <- rpfs %>%
    subset((startoffset < readlen - 5) & (startoffset > 5)) %>%
    mcols(.) %>%
    .[, c("startoffset", "score", "readlen", "phase")] %>%
    as.data.frame() %>% 
    group_by(phase, readlen, startoffset) %>%
    tally(wt = score) %>%
    arrange(startoffset) %>%
    group_by(phase, readlen) %>%
    mutate(gain = cumsum(n)) %>%
    select(p_offset = startoffset, phase, gain)
  # data frame with numbers of reads lost at a specific offset
  lossdf <- rpfs %>%
    subset((.$endoffset < (readlen - 3)) & (.$endoffset > 3)) %>%
    mcols(.) %>%
    .[, c("endoffset", "score", "readlen", "phase")] %>%
    as.data.frame() %>%
    group_by(phase, readlen, endoffset) %>%
    tally(wt = score) %>%
    arrange(endoffset) %>%
    group_by(phase, readlen) %>%
    mutate(loss = cumsum(n)) %>%
    select(p_offset = endoffset, phase, loss) %>%
    mutate(p_offset = p_offset + 3)
  # for a given offset, our gain is the start overlappers whose
  # offset is that or less, minus the stop guys whose offset is
  # less than that
  #
  netdf <- gaindf %>%
    full_join(lossdf, by = c("readlen", "p_offset", "phase")) %>%
    mutate(net = replace_na(gain, 0) - replace_na(loss, 0)) %>%
    group_by(readlen) %>%
    arrange(readlen, p_offset) %>%
    filter(!is.na(net))
  #
  readlens <- rpfs$readlen %>% unique()
  #
  allpos <-
    lapply(readlens, function(rl) {
      expand.grid(
        readlen = rl,
        phase = 0:2,
        p_offset = 2:(rl - 3)
      )
    }) %>%
    bind_rows() %>%
    filter((p_offset %% 3) == 0)
  netdf <- allpos %>% left_join(netdf, 
    by = c("readlen", "p_offset", "phase"))
  netdf <- netdf %>%
    mutate(net = replace_na(net, 0)) %>%
    mutate(twind = net + lag(net) + lead(net))
  netdf <- netdf %>%
    group_by(readlen) %>%
    filter(max(net) > 32)
  best_offsets <- netdf %>%
    filter(!is.na(twind)) %>%
    group_by(readlen) %>%
    dplyr::slice((which.max(twind) - 1):(which.max(twind) + 1))
  #
  uniqueoffsetsfound <- best_offsets %>%
    group_by(readlen, phase) %>%
    filter(n() > 1) %>%
    nrow() %>%
    `==`(0)
  stopifnot(uniqueoffsetsfound)
  #
  best_offsets
}



#' Get Offsets by maximizing p-sites within the CDS
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param input GRanges object with RPFs or a bam file with RPFs
#' @param anno annotation object
#' @param method Method for p-site determination. Defaults to 'cds_incl'
#'
#' @return a dataframe with p-site offsets per readlength
#' @examples
#' cov <- get_readgr(ribobam, anno)
#' offsets_df <- get_offsets(cov, anno)
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#'  data(offsets_df)
#'  testoffsets_df <- get_offsets(rpfs, chr22_anno)

# TODO Add reference to a-p-site paper
#input = ribobam
get_offsets <- function(input, anno, method = "cds_incl") {
  if(is.character(input)&(length(input)==1)){
    message('loading footprint data')
    rpfs <- get_readgr(input, anno)
  }else{
    stopifnot(is(input,'GRanges'))
    rpfs <- input
  }
  #
  if (method == "cds_incl") {
    best_offsets <- get_incl_max_offsets(rpfs, anno)
  } else {
    stop("not yet implemented")
  }
  best_offsets
}
