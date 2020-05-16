# sp3_cols <- c("Protein", "Cleavage_site", "Prob", "Tool")
# sp5_cols <- c("Protein", "Tool", "Classification", "Nterminus", "Cleavage_site", "Prob", "Phase", "X", "Y")
# tmhmm_cols <- c("Protein", "Length", "Exp_AA", "First_60", "Pred_hel", "Topology")

#' Process SignalP-3.1 output
#'
#' Read in a file from SignalP-3.1.
#' Right now, this only handles the hmm setting
#' This function simply reads in the file and adds logical column names.
#' @param file tab-separated output file from SignalP-3.1
#' @keywords signalp
#' @return A tibble
#'
#' @export
read_sp3 <- function(file) {
  cols <- c("Protein", "Cleavage_site", "Prob", "Tool")
  readr::read_tsv(file, col_names = cols)
}

#' Process SignalP-5 output
#'
#' Read in a file from SignalP-5.
#' Right now, it only handles gff output (not the default SignalP-5 setting).
#' This function simply reads in the file and adds logical column names.
#' @param file gff output file from SignalP-5
#' @keywords signalp
#' @return A tibble
#'
#' @export
read_sp5 <- function(file) {
  # columns are from signalp5 gff output
  cols <- c("Protein", "Tool", "Classification", "Nterminus", "Cleavage_site", "Prob", "Phase", "X", "Y")
  readr::read_tsv(file, skip = 1, col_names = cols)
}

#' Process TMHMM-2.0c output
#'
#' Read in a file from TMHMM-2.0c.
#' This function simply reads in the file and adds logical column names
#' @param file gff output file from SignalP-5
#' @keywords tmhmm
#' @return A tibble
#' @seealso \code{/link{pivot_tmhmm}} for conversion to "long" format.
#'
#' @export
read_tmhmm <- function(file, cols=c("Protein", "Length", "Exp_AA", "First_60", "Pred_hel", "Topology")) {
  tmhmm_df <- readr::read_tsv(file, col_names = cols)

  # clean data of text repeated in every row
  tmhmm_clean <- tmhmm_df %>%
    tidyr::separate(Length, c(NA, "Length"), convert = TRUE) %>%
    tidyr::separate(Exp_AA, c(NA, "Exp_AA"), convert = TRUE) %>%
    tidyr::separate(First_60, c(NA, "First_60"), convert = TRUE) %>%
    tidyr::separate(Pred_hel, c(NA, "Pred_hel"), convert = TRUE) %>%
    # if you don't merge the extra pieces of topology, dplyr only keeps the first feat
    tidyr::separate(Topology, c(NA, "Topology"), convert = TRUE, extra = "merge")

  return(tmhmm_clean)
}

#' Tidy TMHMM-2.0c data
#'
#' After reading TMHMM-2.0c output into R, it's still in an unfriendly-format.
#' This tidies the "topology" string in the last column of the output
#' @param x tab-separated TMHMM-2.0c output
#' @keywords tmhmm
#' @return A tibble
#' @seealso \code{/link{read_tmhmm}} for reading input with proper names.
#'
#' @export
#' @examples
#' tmhmm_test <- tibble(Protein = c("A", "B", "C", "D", "E"),
#'                      Info = c(1, 2, 3, 4, 5),
#'                      Topology = c("o", "i", "i4-9",
#'                                   "o10-20i30-50",
#'                                   "i105-205o305-405i505-605"))
#'
#' pivot_tmhmm(tmhmm_test)
#'
#' # This can now be used in ggplot
#' # e.g. How much residue length is TM-helix folds?
#'
#' pivot_tmhmm(tmhmm_test) %>%
#'   mutate(local_in = Stop-Start+1) %>%
#'   group_by(Protein) %>%
#'   mutate(Total_in = sum(local_in, na.rm = TRUE)) %>%
#'   distinct(Protein, .keep_all = TRUE) %>%
#'   select(-Local_in) %>%
#'   ggplot2::ggplot(.) +
#'   ggplot2::geom_histogram(aes(Total_in))
pivot_tmhmm <- function(x) {
  # tidy raw tmhmm output (x) into long format (returns)
  # This was immensely helpful
  # https://www.r-bloggers.com/strsplit-but-keeping-the-delimiter/
  # I thought I needed to use regex for one tmhmm "term",
  # but that didnt end up working anyway
  # topo_regex <- "[io][[:digit:]]+-[[:digit:]]+"
  # solution aws a regex lookahead

  # first get each "term" on its own line
  # then separate each term into columns
  # c(Topology, c("Localization", "Start", "Stop"))
  # Topology is "o" or "i"
  # Start and stop are separated by hyphen

  # I want to split terms into their own line, but avoid the first character
  # so maybe check if i or o AND if previous char was a digit
  tidy_tmhmm <- x %>%
    tidytext::unnest_tokens(Topo_long, Topology,
                            token = stringr::str_split,
                            pattern = "(?=[[:alpha:]])" ) %>%
    dplyr::filter(Topo_long != "") %>%
    tidyr::separate(Topo_long, c("Localization", "pos"), "(?<=[[:alpha:]])") %>%
    tidyr::separate(pos, c("Start", "Stop"), sep = "-",
                    convert = TRUE, fill = "left")

  return(tidy_tmhmm)
}
