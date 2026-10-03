#' Create a block diagonal matrix
#'
#' @param ... Either multiple square matrices, or a single \emph{list} of
#'            square matrices. A single bare matrix is \emph{not} accepted:
#'            with one argument it is coerced via \code{as.list()} into a list
#'            of its scalar elements and rejected. Wrap it in a list instead,
#'            e.g. \code{block_diag(list(m))}.
#'
#' Every input must be a square matrix, and the inputs must either all carry
#' row and column names or all lack them (mixing named and unnamed matrices is
#' an error).
#' @export
block_diag <- function(...) {
    if (...length() > 1) {
        x <- list(...)
    } else {
        x <- as.list(...)
    }
    p <- 0
    rn <- cn <- NULL
    for (m in x) {
        if (!is.matrix(m) || ncol(m) != nrow(m)) {
            stop("Input to `block_diag()` must be a list",
                 "or multiple arguments of square matrices.")
        }
        rn <- c(rn, rownames(m))
        cn <- c(cn, colnames(m))
        p <- p + ncol(m)
    }
    if (length(rn) + length(cn) > 0 && (length(rn) != p || length(cn) != p)) {
        stop("Matrices must all have column and row names, ",
             "or all have no names.")
    }
    out <- matrix(0, nrow = p, ncol = p,
                  dimnames = list(rn, cn))
    last_col <- 0
    for (m in x) {
        ind <- last_col + seq_len(nrow(m))
        out[ind, ind] <- m
        last_col <- last_col + ncol(m)
    }
    out
}