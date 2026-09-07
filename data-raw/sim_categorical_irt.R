# Simulate the categorical-items 2S-PA study -> vignettes/sim_categorical_irt.RDS
#
# Run from the package root:
#   Rscript data-raw/sim_categorical_irt.R           # full (B = 200)
#   Rscript data-raw/sim_categorical_irt.R 20        # quick smoke test (B = 20)
#
# Design: two correlated unit-variance factors F1, F2 with a structural path
#   F2 = gamma* * F1 + zeta,  gamma* = 0.6,
# measured by 5 categorical items each. Items are either binary 2PL or
# polytomous (graded-response, 4 categories). Item thresholds are either
# "balanced" (spread across the latent range) or "skewed" (shifted so most
# responses fall in the lowest category). For each (itemtype, skew, n, rep) we
# estimate the standardized F2 ~ F1 path by three methods:
#   tspa        two-stage path analysis, pooled lavaan stage 2
#   tspa_mx     two-stage path analysis, exact OpenMx stage 2
#   WLSMV       fully joint lavaan with categorical indicators (gold standard)
# and record the standardized estimate, its (delta-method) SE, and convergence.

suppressMessages({
  library(mirt); library(lavaan); library(OpenMx)
  if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE)
  else library(R2spa)
})

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1L) as.integer(args[1L]) else 200L

gamma_star <- 0.6
n_list <- c(50L, 200L, 500L)
a_disc <- 1.5
n_items <- 5L
K <- 4L                       # graded categories (K - 1 = 3 thresholds)
out <- file.path("vignettes", "sim_categorical_irt.RDS")

# ---------------------------------------------------------------- item gen --
# binary 2PL: P(Y = 1 | th) = 1 / (1 + exp(-(a*th - b))); larger b => harder
gen_binary <- function(th, a, b) rbinom(length(th), 1L, 1 / (1 + exp(-(a * th - b))))
# graded response: c_k = P(Y >= k) = 1/(1+exp(-(a*th - d_k))); one uniform draw
# per person, Y = sum_k I(u <= c_k) (categories 0..K-1)
gen_grm <- function(th, a, d) {
  cums <- 1 / (1 + exp(-outer(th, d, function(t, dk) a * t - dk)))
  u <- matrix(runif(nrow(cums)), nrow = nrow(cums), ncol = ncol(cums), byrow = TRUE)
  rowSums(u <= cums)
}

# threshold grids: 5 items, each with 1 (binary) or K-1 (graded) thresholds
# balanced: centered; skewed: shifted high => low marginal responses
thr <- list(
  binary_balanced = seq(-1.5, 1.5, length.out = n_items),
  binary_skewed   = seq(0.5, 3.0, length.out = n_items),
  graded_balanced = sapply(0:(n_items - 1L), function(j) c(-1, 0, 1) + 0.1 * j),
  graded_skewed   = sapply(0:(n_items - 1L), function(j) c(0.5, 1.5, 2.5) + 0.1 * j)
)

gen_items <- function(itemtype, skew, th) {
  key <- paste0(itemtype, "_", skew)
  if (itemtype == "binary") {
    d <- thr[[key]]; sapply(seq_len(n_items), function(j) gen_binary(th, a_disc, d[j]))
  } else {
    D <- thr[[key]]; apply(D, 2, function(dj) gen_grm(th, a_disc, dj))
  }
}

# ------------------------------------------------------- standardized path --
# std.all of the F2 ~ F1 path = gamma * SD(F1) / SD(F2), SD(F2)^2 = v2 + g^2*v1
# where v1 = Var(F1), v2 = residual Var(F2). One formula for all three arms
# (verified to equal lavaan's est.std for the two lavaan arms and, by the same
# algebra, the OpenMx arm). SE by a numeric delta method on (g, v1, v2).
std_path <- function(g, v1, v2) g * sqrt(v1) / sqrt(v2 + g^2 * v1)
se_std_path <- function(p, V, h = 1e-6) {
  f <- function(x) std_path(x[1], x[2], x[3])
  sc <- h * pmax(1, abs(p)); g <- numeric(3)
  for (i in 1:3) { xp <- p; xm <- p; xp[i] <- xp[i] + sc[i]; xm[i] <- xm[i] - sc[i]
    g[i] <- (f(xp) - f(xm)) / (xp[i] - xm[i]) }
  sqrt(max(as.numeric(t(g) %*% V %*% g), 0))
}
subvcov3 <- function(vv, nm) {
  V <- matrix(0, 3, 3, dimnames = list(nm, nm)); idx <- which(nm %in% rownames(vv))
  if (length(idx)) V[idx, idx] <- vv[nm[idx], nm[idx]]; V
}
extract_lavaan <- function(fit) {
  nm <- c("F2~F1", "F1~~F1", "F2~~F2"); cf <- coef(fit)
  if (!all(nm %in% names(cf))) return(c(est = NA_real_, se = NA_real_))
  p <- unname(cf[nm]); c(est = std_path(p[1], p[2], p[3]),
                         se = se_std_path(p, subvcov3(vcov(fit), nm)))
}
extract_mx <- function(fx) {
  nm <- c("m1.A[4,3]", "m1.S[3,3]", "m1.S[4,4]"); cc <- coef(fx)
  if (!all(nm %in% names(cc))) return(c(est = NA_real_, se = NA_real_))
  p <- unname(cc[nm]); c(est = std_path(p[1], p[2], p[3]),
                         se = se_std_path(p, subvcov3(vcov(fx), nm)))
}

seed_for <- function(itemtype, skew, n, b) {
  it <- if (itemtype == "binary") 1L else 2L
  sk <- if (skew == "balanced") 1L else 2L
  as.integer(it * 100000000L + sk * 1000000L + n * 1000L + b)
}

# --------------------------------------------------------------- run one rep --
run_one <- function(itemtype, skew, n, b) {
  set.seed(seed_for(itemtype, skew, n, b))
  e <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, 0, 0, 1 - gamma_star^2), 2))
  th <- cbind(e[, 1], gamma_star * e[, 1] + e[, 2])
  A1 <- gen_items(itemtype, skew, th[, 1]); A2 <- gen_items(itemtype, skew, th[, 2])
  dat <- as.data.frame(cbind(A1, A2)); colnames(dat) <- paste0("i", 1:10)

  rec <- function(method, est, se)
    data.frame(n = n, itemtype = itemtype, skew = skew, method = method,
               est_std = est, se_std = se, converged = !is.na(est),
               ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
               stringsAsFactors = FALSE)

  # stage 1: mirt
  mf <- tryCatch(
    suppressWarnings(mirt(dat, "F1 = 1-5\nF2 = 6-10\nCOV = F1*F2",
                          itemtype = if (itemtype == "binary") "2PL" else "graded")),
    error = function(e) NULL)
  if (is.null(mf)) {
    m <- c("tspa", "tspa_mx", "WLSMV")
    return(do.call(rbind, lapply(m, function(mm) rec(mm, NA_real_, NA_real_))))
  }
  fs <- get_fs(mf)
  est_of <- function(v) c(est = unname(v["est"]), se = unname(v["se"]))

  # stage 2: tspa (pooled)
  r_tspa <- tryCatch({ f <- suppressWarnings(tspa("F2 ~ F1", data = fs));
                       est_of(extract_lavaan(f)) },
                     error = function(e) c(est = NA_real_, se = NA_real_))
  # stage 2: tspa_mx (exact)
  r_mx <- tryCatch({ invisible(capture.output(
                      fx <- suppressWarnings(tspa_mx_model("F2 ~ F1", data = fs))))
                    est_of(extract_mx(fx)) },
                   error = function(e) c(est = NA_real_, se = NA_real_))
  # gold: WLSMV
  ord <- dat; for (j in colnames(ord)) ord[[j]] <- factor(ord[[j]], ordered = TRUE)
  r_w <- tryCatch({ w <- suppressWarnings(sem(
                     "F1 =~ i1+i2+i3+i4+i5\nF2 =~ i6+i7+i8+i9+i10\nF2 ~ F1",
                     data = ord, estimator = "WLSMV"))
                   est_of(extract_lavaan(w)) },
                  error = function(e) c(est = NA_real_, se = NA_real_))

  do.call(rbind, list(
    rec("tspa", unname(r_tspa["est"]), unname(r_tspa["se"])),
    rec("tspa_mx", unname(r_mx["est"]), unname(r_mx["se"])),
    rec("WLSMV", unname(r_w["est"]), unname(r_w["se"]))))
}

# ------------------------------------------------------------------- loop --
cells <- expand.grid(itemtype = c("binary", "graded"),
                     skew = c("balanced", "skewed"),
                     n = n_list, stringsAsFactors = FALSE)
tasks <- expand.grid(cell = seq_len(nrow(cells)), rep = seq_len(B),
                     stringsAsFactors = FALSE)

library(parallel)
nc <- max(1L, min(parallel::detectCores() - 1L, 16L))
run_idx <- function(k) {
  cl <- cells[tasks$cell[k], ]
  run_one(cl$itemtype, cl$skew, cl$n, tasks$rep[k])
}
res_list <- vector("list", nrow(tasks))
t0 <- proc.time()
chunk <- 200L
for (start in seq(1L, nrow(tasks), by = chunk)) {
  idxs <- seq(start, min(start + chunk - 1L, nrow(tasks)))
  res_list[idxs] <- parallel::mclapply(idxs, run_idx, mc.cores = nc, mc.silent = TRUE)
  filled <- which(!vapply(res_list, is.null, logical(1)))
  saveRDS(list(res = do.call(rbind, res_list[filled]), meta = list(B = B)), out)
  message(sprintf("  through task %d/%d  (nc=%d, %.0fs elapsed)",
                  max(idxs), nrow(tasks), nc, (proc.time() - t0)[3]))
}

res <- do.call(rbind, res_list)
# mean marginal item response per condition (a DGP property of the thresholds;
# computed deterministically from a large representative draw)
marg_of <- function(itemtype, skew, n = 5000L) {
  set.seed(99L); th <- rnorm(n)
  mean(rowMeans(as.data.frame(gen_items(itemtype, skew, th))))
}
meta <- list(gamma_star = gamma_star, B = B, n = n_list, seed = 100000L,
             a = a_disc, n_items = n_items, K = K, marginal = c(
               binary_balanced = marg_of("binary", "balanced"),
               binary_skewed = marg_of("binary", "skewed"),
               graded_balanced = marg_of("graded", "balanced"),
               graded_skewed = marg_of("graded", "skewed")))
saveRDS(list(res = res, meta = meta), out)

conv <- aggregate(converged ~ n + itemtype + skew + method, res, mean)
conv$converged <- round(conv$converged, 3)
cat("\nconvergence by condition (fraction of reps):\n")
print(conv)
cat("\nwrote", out, "(", nrow(res), "rows, ",
    (proc.time() - t0)[3], "s )\n")
