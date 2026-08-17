#' Compute interaction indicators for tspa() function
#'
#' This function computes product indicators and the corresponding
#' standard errors according to equation (3) in Hsiao et al. (2021).
#' The double-mean-centering (DMC) strategy is applied in this function,
#' where first-order indicators are mean-centered first and product indicators
#' for latent interaction term(s) are mean-centered again (Lin et al., 2010)
#'
#' @param dat A data frame containing first-order factor score indiactors
#'   with standard error.
#' @param fs_name A vector indicating names of factor scores
#' @param se_fs A vector indicating standard error of factor scores
#' @param loading_fs A vector indicating model-implied loadings of factor scores
#' @param lat_var A vector indicating latent variances
#' @param model An optional string specifying the measurement model
#'              in \code{lavaan} syntax.
#' @return A data frame of product indicators for interaction terms,
#'         with their loadings and standard errors.
#' @references Hsiao, Y.-Y., Kwok, O.-M., & Lai, M. H. C. (2021). Modeling
#'   measurement errors of the exogenous composites from congeneric measures
#'   in interaction models. Structural Equation Modeling: A Multidisciplinary
#'   Journal, 28(2), 250--260. https://doi.org/10.1080/10705511.2020.1782206
#'
#' @importFrom utils combn
#'
#' @export
#'
#' @examples
#' library(lavaan)
#' fs1 <- get_fs(HolzingerSwineford1939,
#'               model = "visual  =~ x1 + x2 + x3",
#'               std.lv = TRUE)
#' fs2 <- get_fs(HolzingerSwineford1939,
#'               model = "textual =~ x4 + x5 + x6",
#'               std.lv = TRUE)
#' fs_dat <- cbind(fs1, fs2)
#' fs_dat2 <- get_fs_int(fs_dat,
#'   fs_name = c("fs_visual", "fs_textual"),
#'   se_fs = c("fs_visual_se", "fs_textual_se"),
#'   loading_fs = c("visual_by_fs_visual",
#'                  "textual_by_fs_textual")
#' )
#' head(fs_dat2[c("fs_visual", "fs_textual", "fs_visual:fs_textual",
#'                "fs_visual:fs_textual_se", "fs_visual:fs_textual_ld")])
get_fs_int <- function(dat, fs_name, se_fs, loading_fs,
                       lat_var = rep_len(1, length.out = length(fs_name)),
                       model = NULL) {
  # Check inputs
  if (!is.data.frame(dat)) {
    stop("'dat' must be a data frame.")
  }
  check_inputs(fs_name, "fs_name", names(dat))
  check_inputs(se_fs, "se_fs", names(dat))
  check_inputs(loading_fs, "loading_fs", names(dat))
  if (!is.numeric(lat_var) || length(lat_var) != length(fs_name)) {
    stop("'lat_var' must be the same length as 'fs_name'.")
  }

  # Create fs pairs
  if (is.null(model)) {
    fs_pairs <- combn(fs_name, 2, simplify = FALSE)
  } else {
    fs_pairs <- lapply(unlist(strsplit(model, split = "\\+")),
      FUN = function(pair) {
        pair_nospace <- trimws(pair)
        unlist(strsplit(pair_nospace, split = ":"))
      }
    )
  }

  # Create PI with SE
  dat_pi <- dat
  for (i in seq_along(fs_pairs)) {
    pair <- fs_pairs[[i]]
    name_i <- paste(pair, collapse = ":")
    # Mean-centered factor product score
    dat_pi[[name_i]] <- dat[[pair[1]]] * dat[[pair[2]]] -
      mean(dat[[pair[1]]] * dat[[pair[2]]])
    # Compute se and loading of the product indicator
    name_se_i <- paste0(name_i, "_se")
    name_ld_i <- paste0(name_i, "_ld")
    se_vars <- se_fs[match(pair, table = fs_name)]
    loading_vars <- loading_fs[match(pair, table = fs_name)]
    lat_vars <- lat_var[match(pair, table = fs_name)]
    dat_pi[[name_se_i]] <- sqrt(
      dat[[loading_vars[1]]]^2 * dat[[se_vars[2]]]^2 * lat_vars[1] +
        dat[[loading_vars[2]]]^2 * dat[[se_vars[1]]]^2 * lat_vars[2] +
        dat[[se_vars[1]]]^2 * dat[[se_vars[2]]]^2
    )
    # Also loadings
    dat_pi[[name_ld_i]] <- dat[[loading_vars[1]]] * dat[[loading_vars[2]]]
  }
  return(dat_pi)
}

# Helper function for input check
check_inputs <- function(input, input_name, names_dat) {
  if (!is.character(input)) {
    stop(paste("'", input_name, "' must be a character vector.", sep = ""))
  }
  if (!all(input %in% names_dat)) {
    missing <- input[!input %in% names_dat]
    stop(paste("The following element(s) in '", input_name,
      "' do(es) not match any column name(s) in 'dat': ",
      paste(missing, collapse = ", "),
      sep = ""
    ))
  }
}
