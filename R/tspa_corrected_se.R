#' First-order correction of sampling covariance for 2S-PA estimates
#' @param tspa_fit A fitted model from [tspa()].
#' @param vfsLT The sampling covariance matrix of `fsL` and `fsT`, which can be
#'              obtained with [get_fs()] with the argument `vfsLT = TRUE`.
#' @param which_free An optional numeric vector indicating which parameters
#'                   in `fsL` and `fsT` are free. The parameters are ordered
#'                   by the `fsL` matrix and the lower-triangular part of
#'                   `fsT`, by columns. For example, for a two-factor model,
#'                   `fsL` and `fsT` are both 2 x 2 matrices, and the
#'                   error covariance between the two factor scores (i.e.,
#'                   the `[2, 1]` element in `fsT`) has an index of 6.
#' @param ... Currently not used.
#' @return A corrected covariance matrix in the same dimension as
#'     `vcov(tspa_fit)`.
#' @export
vcov_corrected <- function(tspa_fit, vfsLT, which_free = NULL, ...) {
    if (is.null(attr(tspa_fit, "fsT"))) {
        stop("corrected vcov requires a tspa model with ",
             "vc and cross-loadings specified.")
    }
    # Direct slot access; avoids lavInspect()'s expensive per-call
    # version check.
    ngrp <- tspa_fit@Data@ngroups
    val_fsL <- attr(tspa_fit, "fsL")
    val_fsT <- attr(tspa_fit, "fsT")
    if (ngrp == 1) {
        if (!is.list(val_fsL)) {
            val_fsL <- list(val_fsL)
        }
        if (!is.list(val_fsT)) {
            val_fsT <- list(val_fsT)
        }
    }
    val_fsT_lower <- lapply(val_fsT, function(x) {
        x[lower.tri(x, diag = TRUE)]
    })
    val_fsLT <- c(unlist(val_fsL), unlist(val_fsT_lower))
    if (is.null(which_free)) {
        which_free <- seq_along(val_fsLT)
    }
    jac <- lavaan::lav_func_jacobian_complex(
        function(x) {
            par <- val_fsLT
            par[which_free] <- x
            counter <- 0
            num_ld <- length(val_fsL[[1]])
            for (g in seq_len(ngrp)) {
                val_fsL[[g]][] <- par[(1:num_ld) + counter]
                counter <- counter + num_ld
            }
            num_ev <- sum(lower.tri(val_fsT[[1]], diag = TRUE))
            for (g in seq_len(ngrp)) {
                val_fsT[[g]][] <- lavaan::lav_matrix_lower2full(
                    par[(1:num_ev) + counter]
                )
                counter <- counter + num_ev
            }
            if (ngrp == 1) {
                val_fsL <- val_fsL[[1]]
                val_fsT <- val_fsT[[1]]
            }
            lavaan::coef(update_tspa(tspa_fit, fsL = val_fsL, fsT = val_fsT))
        },
        x = val_fsLT[which_free],
    )
    vcov(tspa_fit) + jac %*% vfsLT %*% t(jac)
}

update_tspa <- function(tspa_fit, fsL, fsT) {
    call <- attr(tspa_fit, which = "tspa_call")
    if (!missing(fsL)) {
        call$fsL <- eval(call$fsL)
        call$fsL <- fsL
    }
    if (!missing(fsT)) {
        call$fsT <- eval(call$fsT)
        call$fsT <- fsT
    }
    eval(call)
}