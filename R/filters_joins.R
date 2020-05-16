
#' Combine SignalP and TMHMM results
#'
#' Further, filter results on whether or not proteins are secreted.
#' This means right now, TMHMM results must be more comprehensive
#' than SignalP results.
#'
#' @param tmhmm Tibble from read_tmhmm()
#' @param sp Tibble from read_sp3(), sp5 not yet supported.
#' @param is.secreted Also filter the results down to only secretions.
#'
#' @seealso \code{/link{read_tmhmm}}
#'
#' @export
join_sp_to_tmhmm <- function(tmhmm, sp, is.secreted = TRUE) {
  tmhmm_sp <- tmhmm %>%
    pivot_tmhmm() %>%
    left_join(sp, by = "Protein")

  if(is.secreted) {
    tmhmm_sp <- filter(tmhmm_sp, Prob > 0.5)
  }

  return(tmhmm_sp)
}

#' Filter proteins on TMHMM results
#'
#' Removes proteins from joined TMHMM-SignalP data
#'
#' @param x Tibble from internally joined data
#'
#' @sealso \code{/link{join_sp_to_tmhmm}}
#'
#' @export
prune_tm <- function(x) {
  # instead of using filter(any()), checking all TM-helices per prot,
  # just check the stop pos of the most C-terminal helix.
  no_tm_prots <- x %>%
    group_by(Protein) %>%
    mutate(Cterm_TM = max(Stop, na.rm = TRUE)) %>%
    filter(Cterm_TM < 70 | is.na(Cterm_TM)) %>%
    distinct(Protein, .keep_all = TRUE) %>%
    mutate(Cterm_TM = na_if(Cterm_TM, "-Inf"))
}

