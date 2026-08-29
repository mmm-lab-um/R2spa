#' First-order correction of sampling covariance for 2S-PA estimates
#'
#' @details
#' `vcov_corrected()` applies the first-order (delta-method) two-stage
#' approximation: the stage-2 covariance `vcov(tspa_fit)` is augmented by
#' `J %*% vfsLT %*% t(J)`, where `J` is the Jacobian of the stage-2
#' estimates with respect to the selected `fsL`/`fsT` free elements. With the
#' default `engine = "fd"`, `J` is estimated by central differences, one full
#' stage-2 refit on each side of each free element while reusing the base
#' fit's coefficients, so the cost is `2 x (number of free elements)` stage-2
#' refits; with `engine = "analytic"` it is instead evaluated refit-free via
#' a closed form (see the `engine` argument below).
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
#' `"analytic"` is refit-free and deterministic: an influence-function closed
#' form (PLAN 16, sections 2.4 and 4.3), `J = -H^{-1} C`, with `H` (the
#' log-likelihood Hessian over the free params) and `C` (the cross-derivative
#' w.r.t. the fixed `fsL`/`fsT` entries) obtained by central-differencing the
#' analytic log-likelihood score. It covers single- and multi-group models,
#' saturated and restricted (df > 0) structural models, and mean-structure
#' models, and is a pure function of the base fit + `vfsLT` (bit-reproducible,
#' no refits). `"fd"` (central finite differences, one stage-2 refit per side
#' of each free element) is retained as the A/B reference; the analytic path
#' falls back to it only for a shape it cannot handle. The two agree to the
#' finite-difference noise floor whenever `"analytic"` applies.
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
#'              d(thetahat)/d(eta)`. `"analytic"` (the default) uses a
#'              refit-free, deterministic influence-function closed form
#'              (PLAN 16, sections 2.4 and 4.3): `J = -H^{-1} C`, with `H`
#'              (the log-likelihood Hessian over the free params) and `C`
#'              (the cross-derivative w.r.t. the fixed `fsL`/`fsT` entries)
#'              obtained by central-differencing the analytic log-likelihood
#'              score. It covers single- and multi-group models, saturated
#'              and restricted (df > 0) structural models, and mean-structure
#'              models, and is a pure function of the base fit + `vfsLT`
#'              (bit-reproducible, no refits). `"fd"` uses central finite
#'              differences (one stage-2 refit on each side of each free
#'              element) and is retained as the A/B reference; the analytic
#'              path transparently falls back to it only for a shape it cannot
#'              handle (an unrecognised free parameter or unequal per-group
#'              free-param counts). The two agree to the finite-difference
#'              noise floor whenever `"analytic"` applies.
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
                           engine = c("analytic", "fd"), ...) {
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
    if (!is.numeric(which_free) ||
        any(!is.finite(which_free)) ||
        any(which_free != floor(which_free)) ||
        any(which_free < 1) || any(which_free > length(val_fsLT)) ||
        anyDuplicated(which_free)) {
        stop("'which_free' must be a numeric vector of distinct ",
             "whole-number positions in 1:", length(val_fsLT), ".")
    }
    which_free <- as.integer(which_free)
    nfree <- length(which_free)
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
    # Per-group diagonal scale for the PSD perturbation tolerance (the same
    # scaling the base-PSD check above uses), so the padding is not swallowed
    # by floating-point noise when 'fsT' is on a large scale.
    dscale <- vapply(val_fsT, function(x) max(1, max(abs(diag(x)))), 1)
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
            # (Weyl); the extra padding is scaled by the group's diagonal
            # magnitude (as in the base-PSD check), so it is not swallowed
            # by floating-point noise at large scale.
            dmin <- min(eigen(mat$fsT[[g]], symmetric = TRUE,
                              only.values = TRUE)$values)
            if (dmin < -sqrt(qT) * step - 1e-10 * dscale[g]) {
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
    # form is not applicable (multigroup, or a structural model that is not
    # exactly saturated, df > 0). engine = "fd" (the default) is
    # byte-identical to the original central-difference implementation below.
    J <- if (engine == "analytic")
        vcov_jacobian_analytic(tspa_fit, names0, which_free) else NULL
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
# stage-2 model, single- and multi-group (PLAN 16, sections 2.4 and 4.3).
# Returns a (p x nfree) matrix (rows = free structural params in coef order,
# columns = the which_free positions) or NULL when the analytic form is not
# applicable (an unrecognised free parameter, unequal per-group free-param
# counts, or any geometry that fails a guard), in which case the caller falls
# back to the finite-difference Jacobian.
#
# J is the influence function of the stage-2 MLE: by the implicit-function
# theorem, at the MLE  J = -H^{-1} C  where H = d^2 l / dtheta^2 (the
# log-likelihood Hessian over the free params) and C = d^2 l / dtheta deta
# (the cross-derivative w.r.t. the fixed fsL/fsT entries). The multiplier is
# the analytic Hessian, NOT vcov(fit): at an exact MLE the two coincide, but
# vcov(fit) differs from -H^{-1} by O(base-fit residual) and the
# non-saturated C is large, so vcov(fit) C is wrong for restricted models
# (Phase-0 gate: 0.19 vs 0.013 corrected-vcov).
#
# Both H and C are obtained refit-free by central-differencing the analytic
# log-likelihood score s(theta, eta) = d l / d theta (a pure matrix function
# of the implied-covariance geometry; no stage-2 refits, no RNG), so the
# result is deterministic to machine precision. The score covers the
# covariance params (the section 2.4/4.3 form, plus the mean-coupling term
# (n/2)(xbar-mu)' Sinv (dSigma/dtheta) Sinv (xbar-mu) when a mean structure is
# present) and the mean params (n . (Sinv (xbar-mu))_m), so the method handles
# regressions, latent (co)variances, and score means, per group.
vcov_jacobian_analytic <- function(tspa_fit, names0, which_free) {
    ngrp <- tsp_ngroups(tspa_fit)
    if (!isTRUE(tsp_converged(tspa_fit))) return(NULL)
    est <- try(lavaan::lavInspect(tspa_fit, "est"), silent = TRUE)
    if (inherits(est, "try-error")) return(NULL)
    args0 <- attr(tspa_fit, "tspa_args")
    data <- args0$data
    gc <- if (ngrp == 1L) NULL else args0$group
    gns <- if (ngrp == 1L) character(0) else names(est)
    p_total <- length(names0)
    if (p_total %% ngrp != 0L) return(NULL)
    pg <- p_total %/% ngrp
    coef_all <- as.numeric(coef(tspa_fit))
    fsL <- attr(tspa_fit, "fsL"); fsT <- attr(tspa_fit, "fsT")
    if (ngrp == 1L) {
        if (!is.list(fsL)) fsL <- list(fsL)
        if (!is.list(fsT)) fsT <- list(fsT)
    } else if (!is.list(fsL) || !is.list(fsT) ||
        length(fsL) != ngrp || length(fsT) != ngrp) {
        return(NULL)
    }
    q0 <- nrow(fsL[[1L]]); q2 <- q0 * q0; nt <- (q0 * (q0 + 1L)) %/% 2L
    nfull_total <- ngrp * (q2 + nt)
    if (max(which_free) > nfull_total) return(NULL)
    J_full <- matrix(0, p_total, nfull_total)
    for (g in seq_len(ngrp)) {
        e <- if (ngrp == 1L) est else est[[gns[g]]]
        L0 <- e$lambda; psi0 <- e$psi
        q <- nrow(psi0)
        if (q != q0 || nrow(L0) != q || ncol(L0) != q) return(NULL)
        lat <- colnames(L0); score_cols <- rownames(L0)
        if (is.null(lat) || any(is.na(lat))) return(NULL)
        tri <- which(lower.tri(diag(q), diag = TRUE), arr.ind = TRUE)
        nfull <- q * q + nrow(tri)
        sub <- if (ngrp == 1L) data else data[data[[gc]] == gns[g], ]
        x <- try(as.matrix(sub[, score_cols]), silent = TRUE)
        if (inherits(x, "try-error") || any(!is.finite(x))) return(NULL)
        n <- nrow(x); xbar <- colMeans(x)
        S_ml <- crossprod(x - xbar) / n
        idx_g <- (g - 1L) * pg + seq_len(pg)
        th_type <- character(pg); th_i <- integer(pg); th_j <- integer(pg)
        th_score <- integer(pg)
        for (kk in seq_len(pg)) {
            nm <- sub("\\.g[0-9]+$", "", names0[idx_g[kk]])
            if (grepl("~~", nm)) {
                pt <- trimws(strsplit(nm, "~~")[[1]])
                if (length(pt) != 2L) return(NULL)
                th_type[kk] <- "cov"
                th_i[kk] <- match(pt[1], lat); th_j[kk] <- match(pt[2], lat)
            } else if (grepl("~1$", nm)) {
                th_type[kk] <- "mean"
                th_score[kk] <- match(sub("~1$", "", nm), score_cols)
            } else if (grepl("~", nm)) {
                pt <- trimws(strsplit(nm, "~")[[1]])
                if (length(pt) != 2L) return(NULL)
                th_type[kk] <- "reg"
                th_i[kk] <- match(pt[1], lat); th_j[kk] <- match(pt[2], lat)
            } else {
                return(NULL)
            }
            if (th_type[kk] == "mean") {
                if (is.na(th_score[kk])) return(NULL)
            } else if (is.na(th_i[kk]) || is.na(th_j[kk])) {
                return(NULL)
            }
        }
        theta0 <- coef_all[idx_g]
        eta0 <- c(as.vector(fsL[[g]]),
                  as.vector(fsT[[g]][lower.tri(fsT[[g]], diag = TRUE)]))
        # Base implied score means from the partable's observed-variable
        # intercepts (est$nu): FIXED for a single-group model with a mean
        # structure (e.g. fsb supplied), or the free-mean estimates for a
        # multigroup model. Strip the free-mean base values so the central
        # difference re-adds them via tv; without this, a fixed nonzero mean
        # leaks a spurious mean-coupling term into the cov-param score.
        nu0 <- if (!is.null(e$nu) && nrow(e$nu) > 0L)
            as.vector(e$nu[score_cols, 1]) else rep(0, q)
        if (length(nu0) != q) nu0 <- rep(0, q)
        mu_fixed <- nu0
        for (k in which(th_type == "mean"))
            mu_fixed[th_score[k]] <- mu_fixed[th_score[k]] - theta0[k]
        # Analytic score s(theta, eta) for the group (cov + mean params).
        score_grp <- function(tv, ev) {
            beta_of <- function(tv) { b <- matrix(0, q, q)
                for (k in which(th_type == "reg")) b[th_i[k], th_j[k]] <- tv[k]
                b }
            psi_of <- function(tv) { s <- matrix(0, q, q)
                for (k in which(th_type == "cov")) {
                    s[th_i[k], th_j[k]] <- tv[k]; s[th_j[k], th_i[k]] <- tv[k] }
                s }
            Lm <- matrix(ev[seq_len(q * q)], q, q)
            Tm <- matrix(0, q, q)
            for (i in seq_len(nrow(tri)))
                Tm[tri[i, 1], tri[i, 2]] <- ev[q * q + i]
            Tm <- Tm + t(Tm) - diag(diag(Tm))
            mu <- mu_fixed
            for (k in which(th_type == "mean")) mu[th_score[k]] <- tv[k]
            M <- solve(diag(q) - beta_of(tv))
            F <- M %*% psi_of(tv) %*% t(M)
            S <- Lm %*% F %*% t(Lm) + Tm
            Sinv <- solve(S); gS <- Sinv %*% S_ml %*% Sinv - Sinv
            d <- xbar - mu; out <- numeric(length(tv))
            for (k in seq_along(tv)) {
                if (th_type[k] == "mean") {
                    out[k] <- n * (Sinv %*% d)[th_score[k], ]
                } else {
                    if (th_type[k] == "reg")
                        dFk <- M[, th_i[k]] %o% F[th_j[k], ] +
                            F[, th_j[k]] %o% M[, th_i[k]]
                    else
                        dFk <- if (th_i[k] == th_j[k]) M[, th_i[k]] %o% M[, th_i[k]]
                               else M[, th_i[k]] %o% M[, th_j[k]] + M[, th_j[k]] %o% M[, th_i[k]]
                    dSigk <- Lm %*% dFk %*% t(Lm)
                    out[k] <- 0.5 * n * sum(diag(gS %*% dSigk)) +
                        0.5 * n * as.numeric(t(d) %*% Sinv %*% dSigk %*% Sinv %*% d)
                }
            }
            out
        }
        # H = d^2 l / dtheta^2, C = d^2 l / dtheta deta: central differences
        # of the analytic score (refit-free and deterministic).
        h <- 1e-6
        H <- matrix(0, pg, pg)
        for (m in seq_len(pg)) {
            ev <- numeric(pg); ev[m] <- 1
            st <- h * max(1, abs(theta0[m]))
            H[, m] <- (score_grp(theta0 + st * ev, eta0) -
                       score_grp(theta0 - st * ev, eta0)) / (2 * st)
        }
        C <- matrix(0, pg, nfull)
        for (j in seq_len(nfull)) {
            ev <- numeric(nfull); ev[j] <- 1
            st <- h * max(1, abs(eta0[j]))
            C[, j] <- (score_grp(theta0, eta0 + st * ev) -
                       score_grp(theta0, eta0 - st * ev)) / (2 * st)
        }
        Jg <- try(-solve(H, C), silent = TRUE)
        if (inherits(Jg, "try-error")) return(NULL)
        rows <- (g - 1L) * pg + seq_len(pg)
        lcols <- (g - 1L) * q2 + seq_len(q2)
        tcols <- ngrp * q2 + (g - 1L) * nt + seq_len(nt)
        J_full[rows, lcols] <- Jg[, seq_len(q2)]
        J_full[rows, tcols] <- Jg[, (q2 + 1L):nfull]
    }
    J <- J_full[, which_free, drop = FALSE]
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
