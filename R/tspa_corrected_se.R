#' First-order correction of sampling covariance for 2S-PA estimates
#'
#' @details
#' `vcov_corrected()` applies the first-order (delta-method) two-stage
#' approximation: the stage-2 covariance `vcov(tspa_fit)` is augmented by
#' `J %*% vfsLT %*% t(J)`, where `J` is the Jacobian of the stage-2
#' estimates with respect to the selected `fsL`/`fsT` free elements. `J` is
#' estimated by central differences, one full stage-2 refit on each side of
#' each free element while reusing the base fit's coefficients, so the cost
#' is `2 x (number of free elements)` stage-2 refits.
#'
#' The correction is **partial by design**. It propagates only the *sampling*
#' uncertainty of the stage-1 estimates of `fsL` (loadings/cross-loadings)
#' and `fsT` (score error variance-covariance, including its off-diagonal).
#' It does **not** account for the sampling uncertainty of the factor-score
#' *values* (the `fs_<name>` columns fitted in stage 2), their standard
#' errors (`se_fs`), or the score intercepts (`fsb`); these are held fixed.
#' The two-stage approximation also treats the stages as independent, so the
#' cross-covariance between the stage-1 estimates and the stage-2 data is not
#' modelled.
#'
#' The same correction can also be requested in place via
#' `tspa(..., corrected_se = TRUE, vfsLT = <matrix>)`; `tspa()` then sets the
#' `tspa_corrected = TRUE` attribute on the returned fit, and
#' `standardizedSolution()` on that fit reports corrected standard SEs. Such
#' an already-corrected fit is rejected by `vcov_corrected()` (the
#' correction is never applied twice).
#'
#' The `engine` argument selects how the Jacobian is evaluated. The default
#' `"fd"` uses central finite differences (one stage-2 refit on each side of
#' each free element). `"analytic"` instead uses a refit-free,
#' deterministic closed form for the (approximately) saturated single-group
#' case (PLAN 16, section 2.4) and transparently falls back to `"fd"`
#' otherwise (multigroup, or a restricted structural model that does not
#' approximately saturate the score covariance); the two agree to the
#' finite-difference noise floor whenever `"analytic"` applies, and the
#' analytic result is bit-reproducible.
#'
#' @param tspa_fit A fit from [tspa()] with `fsT` and `fsL` supplied
#'              (multi-factor measurement model, single- or multi-group),
#'              so that it carries the `fsT`, `fsL`, and `tspa_args`
#'              attributes. A fit corrected in place via
#'              `tspa(corrected_se = TRUE)` (attribute
#'              `tspa_corrected = TRUE`) is rejected.
#' @param vfsLT The sampling covariance matrix of the free `fsL`/`fsT`
#'              elements, taken (or sub-matrixed) from the `vfsLT` attribute
#'              of a [get_fs()] result fitted with `vfsLT = TRUE`. Its row and
#'              column order is the same `fsL`-then-`fsT` order described by
#'              `which_free`. When `which_free = NULL` (the default) this is
#'              the full square matrix returned by `get_fs()`; when
#'              `which_free` is a subset of length `k`, this **must** be the
#'              matching `k x k` principal submatrix, i.e.
#'              `vfsLT_full[which_free, which_free]`.
#' @param which_free An optional numeric vector of positions selecting which
#'                   `fsL`/`fsT` elements to treat as free (and therefore to
#'                   propagate through the Jacobian). Positions run over the
#'                   `fsL` matrix (column-major order) followed by the
#'                   lower-triangular part of `fsT` including the diagonal
#'                   (column-major order). For a two-factor model `fsL` and
#'                   `fsT` are both 2 x 2: `fsL` occupies positions 1:4 and
#'                   `fsT` positions 5:7, so the two error variances are 5
#'                   and 7 and the error covariance between the two factor
#'                   scores (the `[2, 1]` element of `fsT`) is position 6.
#'                   For a multigroup fit the positions run per group, in
#'                   group order: group 1's full `fsL` (column-major), then
#'                   group 2's full `fsL`, ..., then group 1's
#'                   lower-triangular `fsT` (column-major), then group 2's,
#'                   i.e. all loadings across groups first, then all
#'                   error-variance elements (matching the order of the
#'                   `vfsLT` attribute from `get_fs(vfsLT = TRUE)`).
#'                   `NULL` (the default) treats every per-group `fsL`/`fsT`
#'                   element as free. A non-`NULL` `which_free` of length
#'                   `k` requires `vfsLT` to be the matching `k x k`
#'                   principal submatrix (see `vfsLT`).
#' @param engine The engine used to evaluate the Jacobian `J =
#'              d(thetahat)/d(eta)`. `"fd"` (the default) uses central
#'              finite differences (one stage-2 refit on each side of each
#'              free element). `"analytic"` uses a refit-free, deterministic
#'              closed form for the (approximately) saturated single-group
#'              case (PLAN 16, section 2.4) and transparently falls back to
#'              `"fd"` otherwise (multigroup, or a structural model that does
#'              not approximately saturate the score covariance). The two
#'              engines agree to the finite-difference noise floor whenever
#'              `"analytic"` applies.
#' @param ... Currently not used.
#' @return A corrected covariance matrix in the same dimension as
#'     `vcov(tspa_fit)` (symmetric).
#' @examples
#' library(lavaan)
#'
#' # Two-factor model, Bartlett scoring. The fsL/fsT free elements run over
#' # positions 1:7 (fsL -> 1:4, fsT -> 5:7); the two error variances are 5
#' # and 7, the score-error covariance is 6.
#' fs <- get_fs(PoliticalDemocracy,
#'              model = "ind60 =~ x1 + x2 + x3
#'                       dem60 =~ y1 + y2 + y3 + y4",
#'              method = "Bartlett", vfsLT = TRUE)
#' fit <- tspa("dem60 ~ ind60", data = fs,
#'             fsT = attr(fs, "fsT"), fsL = attr(fs, "fsL"))
#' vfsLT <- attr(fs, "vfsLT")
#'
#' # Propagate only the two error variances (positions 5, 7): with a
#' # non-NULL which_free, pass the matching principal submatrix of vfsLT.
#' vc <- vcov_corrected(fit, vfsLT = vfsLT[c(5, 7), c(5, 7)],
#'                      which_free = c(5, 7))
#' sqrt(diag(vc))   # corrected standard errors
#' @export
vcov_corrected <- function(tspa_fit, vfsLT, which_free = NULL,
                           engine = c("fd", "analytic"), ...) {
    engine <- match.arg(engine)
    if (is.null(attr(tspa_fit, "fsT"))) {
        stop("corrected vcov requires a tspa() fit with ",
             "'fsT' and 'fsL' supplied (multi-factor measurement model).")
    }
    if (isTRUE(attr(tspa_fit, "tspa_corrected"))) {
        stop("the fit is already SE-corrected (tspa_corrected = TRUE; from ",
             "tspa(corrected_se = TRUE)); pass the uncorrected fit so the ",
             "correction is not applied twice.")
    }
    args0 <- attr(tspa_fit, "tspa_args")
    if (is.null(args0)) {
        stop("the 'tspa_args' attribute is missing from the fit; ",
             "refit with the current R2spa version.")
    }
    # Consolidated through the lavaan compat boundary (single point of
    # coupling; slot access with lavInspect fallback).
    ngrp <- tsp_ngroups(tspa_fit)
    val_fsL <- attr(tspa_fit, "fsL")
    val_fsT <- attr(tspa_fit, "fsT")
    if (ngrp == 1) {
        val_fsL <- if (is.list(val_fsL)) val_fsL else list(val_fsL)
        val_fsT <- if (is.list(val_fsT)) val_fsT else list(val_fsT)
    } else {
        if (!is.list(val_fsL) || !is.list(val_fsT) ||
            length(val_fsL) != ngrp || length(val_fsT) != ngrp) {
            stop("multigroup fits require 'fsL' and 'fsT' attribute ",
                 "lists of length ", ngrp, ".")
        }
    }
    for (g in seq_len(ngrp)) {
        if (!is.matrix(val_fsL[[g]]) || !is.matrix(val_fsT[[g]]) ||
            !identical(dim(val_fsL[[g]]), dim(val_fsL[[1]])) ||
            !identical(dim(val_fsT[[g]]), dim(val_fsT[[1]])) ||
            !identical(dimnames(val_fsT[[g]]), dimnames(val_fsT[[1]]))) {
            stop("per-group 'fsL'/'fsT' dimensions and dimnames must be ",
                 "identical; group ", g, " differs from group 1.")
        }
        if (!all(is.finite(val_fsL[[g]])) || !all(is.finite(val_fsT[[g]]))) {
            stop("'fsL'/'fsT' must contain only finite values ",
                 "(group ", g, ").")
        }
    }
    # Parameter layout: [fsL_g1 (full, column-major), ..., fsL_gn,
    #                   fsT_g1 (lower triangle, column-major), ..., fsT_gn]
    # — the same order as the 'vfsLT' attribute produced by get_fs().
    val_fsLT <- c(unlist(val_fsL), unlist(lapply(
        val_fsT, function(x) x[lower.tri(x, diag = TRUE)])))
    if (is.null(which_free)) {
        which_free <- seq_along(val_fsLT)
    }
    which_free <- as.integer(which_free)
    nfree <- length(which_free)
    if (anyNA(which_free) || anyDuplicated(which_free) ||
        any(which_free < 1L) || any(which_free > length(val_fsLT))) {
        stop("'which_free' must contain distinct positions in ",
             "1:", length(val_fsLT), ".")
    }
    # Fail-fast input guards (before any refit is spent):
    if (!is.matrix(vfsLT) || nrow(vfsLT) != ncol(vfsLT) ||
        nrow(vfsLT) != nfree) {
        stop("'vfsLT' must be a (", nfree, " x ", nfree, ") matrix ",
             "matching 'which_free'; when 'which_free' is a subset, pass ",
             "the corresponding principal submatrix of the full 'vfsLT' ",
             "attribute, i.e. 'vfsLT[which_free, which_free]'.")
    }
    if (!isTRUE(all.equal(vfsLT, t(vfsLT)))) {
        stop("'vfsLT' must be symmetric.")
    }
    if (!all(is.finite(vfsLT))) {
        stop("'vfsLT' must contain only finite values.")
    }
    # The base fsT must be positive semi-definite (a get_fs() attribute
    # always is; this catches hand-built or corrupted inputs early).
    for (g in seq_len(ngrp)) {
        emin <- min(eigen(val_fsT[[g]], symmetric = TRUE,
                          only.values = TRUE)$values)
        if (emin < -1e-8 * max(1, max(abs(diag(val_fsT[[g]]))))) {
            stop("'fsT' is not positive semi-definite (group ", g,
                 ", smallest eigenvalue ", signif(emin, 3), ").")
        }
    }
    num_ld <- length(val_fsL[[1]])
    num_ev <- sum(lower.tri(val_fsT[[1]], diag = TRUE))
    qT <- nrow(val_fsT[[1]])
    # theta vector -> (fsL, fsT) per group, preserving dimnames.
    assemble <- function(par) {
        counter <- 0L
        Ls <- val_fsL
        for (g in seq_len(ngrp)) {
            Ls[[g]][, ] <- par[(1:num_ld) + counter]
            counter <- counter + num_ld
        }
        Ts <- vector("list", ngrp)
        for (g in seq_len(ngrp)) {
            Ts[[g]] <- tsp_tri2full_colmajor(
                par[(1:num_ev) + counter], qT, dimnames(val_fsT[[g]]))
            counter <- counter + num_ev
        }
        list(fsL = Ls, fsT = Ts)
    }
    # The base reconstruction must round-trip exactly (guards the
    # vectorization wiring; both sides are column-major by construction).
    mat0 <- assemble(val_fsLT)
    if (!all(mapply(function(L, L0) isTRUE(all.equal(L, L0)),
                    mat0$fsL, val_fsL)) ||
        !all(mapply(function(T, T0) isTRUE(all.equal(T, T0)),
                    mat0$fsT, val_fsT))) {
        stop("internal error: fsL/fsT reconstruction round-trip failed.")
    }
    coef0 <- coef(tspa_fit)
    if (length(coef0) == 0L) {
        stop("the base fit has no free parameters; nothing to correct.")
    }
    names0 <- names(coef0)
    # Refit the stage-2 model at a perturbed theta via the self-contained
    # recorded arguments (no environment dependence); only the estimates
    # are needed, so stage-2 SEs are skipped.
    refit_coef <- function(par, k, step) {
        mat <- assemble(par)
        for (g in seq_len(ngrp)) {
            # A PSD base can drop by at most sqrt(q) * step per eigenvalue
            # (Weyl); anything beyond that means the assembly is broken.
            dmin <- min(eigen(mat$fsT[[g]], symmetric = TRUE,
                              only.values = TRUE)$values)
            if (dmin < -sqrt(qT) * step - 1e-10) {
                stop("the perturbed 'fsT' (group ", g, ", theta index ",
                     k, ") is not positive semi-definite; the ",
                     "first-order correction is undefined here. Try ",
                     "restricting 'which_free' or check the stage-1 fit.")
            }
        }
        a <- args0
        if (ngrp == 1) {
            a$fsL <- mat$fsL[[1]]
            a$fsT <- mat$fsT[[1]]
        } else {
            a$fsL <- mat$fsL
            a$fsT <- mat$fsT
        }
        a$se <- "none"
        # Refits must be plain stage-2 fits (no re-correction).
        a$corrected_se <- FALSE
        fit_i <- do.call(tspa, a)
        check_refit_convergence(fit_i, k)
        c_i <- coef(fit_i)
        if (!identical(names(c_i), names0)) {
            stop("internal error: the stage-2 refit changed the ",
                 "free-parameter set (theta index ", k, ").")
        }
        c_i
    }
    # Jacobian J = d(thetahat)/d(eta). engine = "analytic" evaluates it
    # refit-free via the saturated closed form (PLAN 16, section 2.4) and
    # silently falls back to the finite-difference engine when the analytic
    # form is not applicable (multigroup, or a non-saturated structural
    # model). engine = "fd" (the default) is byte-identical to the original
    # central-difference implementation below.
    J <- if (engine == "analytic")
        vcov_jacobian_analytic(tspa_fit, args0, names0, which_free) else NULL
    if (is.null(J)) {
        # Central-difference step: verified stable over h in 1e-5..1e-7 for
        # the package's stage-2 models (the stage-2 MLE is strongly curved,
        # so the h = 1e-4 truncation error is already ~0.3%; optimizer jitter
        # at h = 1e-5 is ~1e-8 relative and negligible).
        h0 <- 1e-5
        J <- matrix(0, nrow = length(coef0), ncol = nfree)
        rownames(J) <- names0
        colnames(J) <- as.character(which_free)
        for (k in seq_len(nfree)) {
            p <- which_free[k]
            step <- h0 * max(1, abs(val_fsLT[p]))
            e <- numeric(length(val_fsLT))
            e[p] <- 1
            c_plus <- refit_coef(val_fsLT + step * e, k, step)
            c_minus <- refit_coef(val_fsLT - step * e, k, step)
            J[, k] <- (c_plus - c_minus) / (2 * step)
        }
    }
    vcov(tspa_fit) + J %*% vfsLT %*% t(J)
}

# Convergence gate for theta-perturbation refits: a non-converged refit
# means the first-order correction is undefined at that perturbation, so
# stop loudly (naming the theta index) rather than feed a garbage column
# into the Jacobian.
check_refit_convergence <- function(fit, k) {
    if (!isTRUE(tsp_converged(fit))) {
        stop("stage-2 refit for Jacobian element ", k,
             " did not converge; the first-order correction is ",
             "undefined at this perturbation. Check that 'fsT' stays ",
             "well-conditioned and positive definite, or restrict ",
             "'which_free'.")
    }
    invisible(TRUE)
}

# Analytic influence-function Jacobian J = d(thetahat)/d(eta) for the
# saturated single-group stage-2 model (PLAN 16, section 2.4). Returns a
# p x nfree matrix (rows = free structural params in coef order, columns =
# the which_free positions) or NULL when the analytic form is not applicable
# (multigroup, a non-saturated structural model, or any geometry the closed
# form does not cover), in which case the caller falls back to the
# finite-difference Jacobian.
#
# The saturated cross-Hessian is Hessian-free (no log-likelihood Hessian, no
# finite differences, no refits):
#   H_theta_eta[k, j] = -(n/2) * tr( Sinv (dSigma/deta_j) Sinv (dSigma/dtheta_k) )
#   J = vcov(fit) %*% H_theta_eta
# where Sigma = L F L' + T, F = (I-beta)^{-1} psi (I-beta)^{-1}', Sinv is
# the inverse of Sigma, and the dSigma/dtheta_k, dSigma/deta_j are the
# analytic derivatives below. It is exact when the structural part saturates
# the score covariance (Sigma = the sample score covariance at the MLE), which
# is checked against the data before the closed form is used.
vcov_jacobian_analytic <- function(tspa_fit, args0, names0, which_free) {
    if (tsp_ngroups(tspa_fit) != 1L) return(NULL)
    est <- try(lavaan::lavInspect(tspa_fit, "est"), silent = TRUE)
    if (inherits(est, "try-error")) return(NULL)
    L <- est$lambda; psi <- est$psi; beta <- est$beta; T0 <- est$theta
    q <- nrow(psi)
    if (nrow(L) != q || ncol(L) != q) return(NULL)
    lat <- rownames(psi)
    if (is.null(lat) || any(is.na(lat))) return(NULL)
    # Observation count (the closed form scales by n).
    data0 <- args0$data
    if (is.null(data0)) return(NULL)
    n <- nrow(data0)
    # Sample covariance of the score columns, taken from lavaan's own
    # observed covariance (the ML estimator, consistent with the fit's
    # objective — the n vs n-1 scaling, weights, and FIML handled exactly as
    # in the fit), used only for the near-saturation gate below. Using
    # stats::cov() here would inject a systematic n-1 scaling bias into the
    # gate (a saturated fit then reads as non-saturated at ~ (1/n) * scale).
    cov_ov <- try(lavaan::lavInspect(tspa_fit, "cov.ov"), silent = TRUE)
    if (inherits(cov_ov, "try-error") ||
        !all(rownames(L) %in% rownames(cov_ov))) return(NULL)
    S_ml <- cov_ov[rownames(L), rownames(L), drop = FALSE]
    if (!all(is.finite(S_ml))) return(NULL)
    # Geometry: F = (I-beta)^{-1} psi (I-beta)^{-1}' (full latent cov, not
    # est$psi, which is exogenous-only), Sigma = L F L' + T.
    M <- try(solve(diag(q) - beta), silent = TRUE)
    if (inherits(M, "try-error")) return(NULL)
    F <- M %*% psi %*% t(M)
    Sigma <- L %*% F %*% t(L) + T0
    Sinv <- try(solve(Sigma), silent = TRUE)
    if (inherits(Sinv, "try-error")) return(NULL)
    # Near-saturation gate. The closed form is exact when the implied score
    # covariance Sigma equals the sample score covariance S_ml at the MLE
    # (the structural part saturates the score covariance); it degrades
    # gracefully when the fit is only close. A restricted structural model
    # that cannot reproduce S_ml leaves the second-order term in the
    # cross-Hessian non-negligible, so fall back to the FD there. Because
    # S_ml is lavaan's own ML observed covariance, a (near-)saturated fit
    # reads as ~0 while a clearly restricted one reads as ~0.5, so the
    # relative threshold 0.1 * mean(diag(S_ml)) cleanly separates the two.
    if (max(abs(Sigma - S_ml)) > 0.1 * mean(diag(S_ml))) return(NULL)
    # Delta_theta: one q x q matrix per free structural param (coef order).
    dFb <- function(i, j) M[, i] %o% F[j, ] + F[, j] %o% M[, i]
    dFp <- function(k) M[, k] %o% M[, k]
    p <- length(names0)
    dSigTheta <- vector("list", p)
    for (kk in seq_len(p)) {
        nm <- names0[kk]
        if (grepl("~~", nm)) {
            # latent (co)variance: "var~~var". Only the diagonal (a free
            # residual/exogenous variance) is covered; an off-diagonal free
            # psi covariance is out of scope for the closed form.
            parts <- trimws(strsplit(nm, "~~")[[1]])
            if (length(parts) != 2L || parts[1] != parts[2]) return(NULL)
            k2 <- match(parts[1], lat)
            if (is.na(k2)) return(NULL)
            dSigTheta[[kk]] <- L %*% dFp(k2) %*% t(L)
        } else if (grepl("~", nm)) {
            # structural regression "lhs~rhs"
            parts <- trimws(strsplit(nm, "~")[[1]])
            if (length(parts) != 2L) return(NULL)
            i2 <- match(parts[1], lat); j2 <- match(parts[2], lat)
            if (is.na(i2) || is.na(j2)) return(NULL)
            dSigTheta[[kk]] <- L %*% dFb(i2, j2) %*% t(L)
        } else {
            return(NULL)  # an unrecognised free parameter: use the FD
        }
    }
    # Delta_eta: q^2 loading + lower-tri error entries, in the fsL-then-fsT
    # (column-major) order of the vfsLT / which_free layout.
    FL <- F %*% t(L); LF <- L %*% F
    dSig_L <- function(r, c2) { ei <- matrix(0, q, 1); ei[r, 1] <- 1
        ei %*% FL[c2, ] + LF[, c2] %*% t(ei) }
    dSig_T <- function(r, c2) { E <- matrix(0, q, q); E[r, c2] <- 1; E[c2, r] <- 1; E }
    tri <- which(lower.tri(diag(q), diag = TRUE), arr.ind = TRUE)
    nfull <- q * q + nrow(tri)
    if (max(which_free) > nfull) return(NULL)
    dSigEta <- vector("list", nfull)
    idx <- 0L
    for (c2 in seq_len(q)) for (r in seq_len(q)) {
        idx <- idx + 1L; dSigEta[[idx]] <- dSig_L(r, c2)
    }
    for (t in seq_len(nrow(tri))) {
        idx <- idx + 1L; dSigEta[[idx]] <- dSig_T(tri[t, 1], tri[t, 2])
    }
    # Cross-Hessian, then J = V %*% H_theta_eta (subset to which_free).
    Hte <- matrix(0, nrow = p, ncol = nfull)
    for (kk in seq_len(p)) for (jj in seq_len(nfull))
        Hte[kk, jj] <- -0.5 * n * sum(diag(Sinv %*% dSigEta[[jj]] %*% Sinv %*%
                                        dSigTheta[[kk]]))
    J <- vcov(tspa_fit) %*% Hte[, which_free, drop = FALSE]
    rownames(J) <- names0
    colnames(J) <- as.character(which_free)
    J
}

# Inverse of m[lower.tri(m, diag = TRUE)]: fills the lower triangle (incl.
# diagonal) column-major and mirrors it to the upper triangle, matching the
# producer-side convention used throughout the package (and in the 'vfsLT'
# attribute of get_fs()). lavaan::lav_mat_lower2full() reads the vector in
# ROW-major order, so it is deliberately not used here.
tsp_tri2full_colmajor <- function(x, q, dn = NULL) {
    m <- q * (q + 1L) / 2L
    if (length(x) != m) {
        stop("expected ", m, " values for the ", q, " x ", q,
             " lower triangle; got ", length(x), ".")
    }
    out <- matrix(0, q, q)
    out[lower.tri(out, diag = TRUE)] <- x
    out <- out + t(out) - diag(diag(out))
    dimnames(out) <- dn
    out
}
