#' Two-Stage Path Analysis
#'
#' Fit a two-stage path analysis (2S-PA) model.
#'
#' @details
#' When the factor-score attributes are heterogeneous across the units of a
#' group, `tspa()` first re-expresses them as long-form, individual-specific
#' values and then reduces them to a single representative set per group;
#' that pooled set is what feeds the stage-2 model and is attached to the
#' returned fit. The heterogeneous cases are per-pattern values from a group
#' fitted with missing data (`missing = "fiml"`), where `fsL`/`fsT`/`fsb` are
#' per-group lists of one matrix/vector per observed-indicator pattern, and
#' per-cluster values from a `merMod` fit, where `fsL`/`fsT` are 3-D arrays
#' (one slice per cluster). Pooling, rather than fitting each pattern as its
#' own tiny stage-2 sub-group, keeps small (possibly near-empty) patterns from
#' making the measurement model under-identified or numerically fragile. The
#' default `reduce` is `"mean"`, a convex combination of positive
#' semi-definite matrices, so the pooled `fsT` stays positive semi-definite;
#' the opt-in `"median"` trades that guarantee for robustness and emits a
#' warning when the pooled `fsT` is not positive semi-definite. For
#' homogeneous inputs the reduction is a no-op, so complete-data behavior is
#' unchanged.
#'
#' ## Derivation from a [get_fs()] result
#'
#' `tspa()` silently derives its measurement inputs from a [get_fs()]
#' result passed as `data`, so the canonical call is
#' `tspa(model, data = get_fs(...))`. Two derivations exist, in this
#' order, and explicit arguments always win over them:
#'
#' - *Multi-factor* — fires when `fsT`, `fsL`, and `se_fs` are all
#'   omitted and `data` carries both `fsT` and `fsL` attributes; those
#'   attributes (and the `fsb` attribute when present) become the
#'   measurement inputs. Derivation is provenance-gated: `data` must
#'   actually resolve as a [get_fs()] result (it carries the
#'   `fs_pattern` attribute, or a `merMod` 3-D / `mirt` per-observation
#'   shape), so a hand-rolled data frame with plain matrix `fsT`/`fsL`
#'   attributes but no such provenance is not derived and is rejected
#'   with an informative error. When both forms are available and
#'   nothing is passed, this (attribute) form wins.
#' - *Single-factor* — fires when `fsT` is still `NULL` after that and
#'   `se_fs` was omitted; `se_fs` is then built from the data's own
#'   `fs_<v>` score columns that carry a matching numeric `fs_<v>_se`
#'   column (simple names only; product-score columns
#'   `fs_<v1>:<v2>` are not derived). This is the path taken by a
#'   `cbind()`'d [get_fs()] result, because `cbind()` drops the
#'   attributes and leaves only the score/SE columns.
#'
#' A supplied `se_fs` (even an empty `list()`) suppresses the
#' multi-factor derivation, keeping the single-factor path.
#'
#' In the single-factor derivation with a group column, the derived
#' `se_fs` carries one row per group, in the order the group column's
#' values first appear in the data (lavaan's `group =` order for both
#' character and factor columns); a within-group-constant SE column
#' keeps its constant and a varying column (FIML missing data) is
#' reduced by `reduce`, idempotent with the FIML pooling that may run on
#' the result. The group column is the data's `group_col` attribute (if
#' that column exists), else the `group =` argument (if it names a
#' column of the data), else a literal `group` column. A `cbind()`'d
#' multi-group frame with no such group signal (no `group_col` attribute
#' — `cbind()` drops it — and no `group =` argument) is silently fitted
#' as single-group with row-mean SEs, so pass `group =` to control the
#' stage-2 grouping.
#'
#' The multi-factor derivation also picks up the `fsb` (intercept)
#' attribute, so a derived fit includes the stage-2 intercept
#' constraints. On `std.lv = TRUE` data the derived intercepts are zero
#' (the default `regression`/`Bartlett` scores), so the derived fit's
#' free estimates and standard errors equal the corresponding values of
#' the explicit no-`fsb` form (the intercept rows are simply fixed); the
#' differences are that its `tspaModel` string carries an extra intercept
#' block and that, for multigroup fits (where lavaan enables a mean
#' structure), the no-`fsb` form additionally estimates the factor-score
#' intercepts freely. On data with nonzero `fsb` (e.g. a `mirt` fit) even
#' the free estimates differ from the no-`fsb` form.
#'
#' ## Product-score auto-compute (`product = TRUE`)
#'
#' With `product = TRUE`, a model latent that names the product of two of
#' the model's factor scores is treated as a latent interaction measured
#' by the double-mean-centered product indicator of the two scores. The
#' pair may be named by concatenation (`xm` for `x` and `m`) or in
#' lavaan's interaction syntax (`x:m`); the interaction-syntax form is
#' rewritten to the concatenated render name, because in the generated
#' model `x:m` would be parsed by lavaan as an interaction of the (latent)
#' variables. An `a:b` token whose parts are not both factor scores (e.g.
#' `x:g` with `g` an observed covariate) is not claimed and is passed
#' through to lavaan as an ordinary interaction. The product columns are
#' computed on the fly when absent ([compute_fs_prod()] from the data's
#' stage-1 attributes; pre-existing `fs_a:fs_b` columns — in either
#' orientation — are used as-is, so a [get_fs()] result with `product` set
#' works too) and the product is wired into the stage-2 measurement model:
#' in the single-factor path the product SE joins `se_fs` (per-group
#' pooled by `reduce`, the same convention as the score SEs; an explicit
#' product SE may be keyed by either the render name `xm` or the token
#' `x:m`), and in the multi-factor path the product latent gets a fixed
#' loading `gamma` and fixed error variance `se_P^2` from the (pooled)
#' `fsL`/`fsT`/`psi`. Rejected with an informative error (v1): multigroup
#' models, `corrected_se = TRUE`, and data without stage-1 attributes that
#' lacks the product columns (e.g. a `cbind()`'d [get_fs()] result —
#' compute the product columns up front or pass the un-`cbind()`'d
#' result). A model variable matching two different factor-score pairs,
#' the same pair named twice (`x:m` and `xm`), or a render name colliding
#' with another model variable, is an error.
#'
#' @param model A string variable describing the structural path model,
#'              in \code{lavaan} syntax.
#' @param data A data frame containing factor scores. When `data` is a
#'              [get_fs()] result and the corresponding arguments are
#'              omitted, the measurement inputs are derived from it: the
#'              `fsT`/`fsL`/`fsb` attributes for a multi-factor fit, or,
#'              for a `cbind()`'d frame whose attributes `cbind()` drops,
#'              the `fs_<v>`/`fs_<v>_se` score columns for a
#'              single-factor fit (see `Details`).
#' @param reliability A numeric vector representing the reliability indexes
#'                    of each latent factor. Currently \code{tspa()} does not
#'                    support the reliability argument. Please use \code{se}.
#' @param se Deprecated to avoid conflict with the argument of the same name
#'           in [lavaan::lavaan()].
#' @param se_fs A numeric vector representing the standard errors of each
#'              factor score variable for single-group 2S-PA. A list or data
#'              frame storing the standard errors of each group in each latent
#'              factor for multigroup 2S-PA. An explicit `se_fs` always wins
#'              over derivation from a [get_fs()] result and, even when
#'              empty (`list()`), suppresses the multi-factor (attribute)
#'              derivation; the single-factor `se_fs` is derived from the
#'              data's `fs_<v>_se` columns only when this argument is
#'              omitted.
#' @param fsT An error variance-covariance matrix of the factor scores, which
#'            can be obtained from the output of [get_fs()] using `attr()`
#'            with the argument `which = "fsT"`. When a group was fitted with
#'            missing data (`missing = "fiml"`), the attribute carries
#'            per-pattern values (a per-group list of one matrix per
#'            observed-indicator pattern); for a `merMod` fit it is a 3-D
#'            per-cluster array; for a `mirt` fit it is a per-observation list
#'            (one matrix per row, marked `mirt_per_obs`). Values of these
#'            per-unit shapes are reduced to a single representative per-group
#'            matrix by `reduce`; the pooled (not the nested/per-cluster)
#'            matrix is attached to the returned fit as the `fsT` attribute.
#'            When omitted, `fsT` is derived from the `fsT` attribute of a
#'            [get_fs()] result passed as `data` (see `Details`); an
#'            explicit `fsT` always wins.
#' @param fsL A matrix of loadings and cross-loadings from the
#'            latent variables to the factor scores `fs`, which
#'            can be obtained from the output of [get_fs()] using
#'            `attr()` with the argument `which = "fsL"`.
#'            For details see the Multi-Factor Measurement Model vignette:
#'            `vignette("Multi-Factor Measurement Model", package = "R2spa")`.
#'            As with `fsT`, per-pattern (FIML missing data), per-cluster
#'            (merMod), and per-observation (mirt) values are supported and
#'            reduced per group by `reduce`; the pooled (not the nested)
#'            matrix is attached to the returned fit as the `fsL` attribute.
#'            When omitted, `fsL` is derived from the `fsL` attribute of a
#'            [get_fs()] result passed as `data` (see `Details`); an
#'            explicit `fsL` always wins.
#' @param fsb A vector of intercepts for the factor scores `fs`, which can
#'            be obtained from the output of [get_fs()] using `attr()`
#'            with the argument `which = "fsb"`. As with `fsT`, per-pattern
#'            (FIML missing data), per-cluster (merMod), and per-observation
#'            (mirt) values are supported
#'            and reduced per group by `reduce`; the pooled (not the
#'            nested/per-cluster) vector is used for the stage-2 intercept
#'            constraints. When omitted, `fsb` is derived from the `fsb`
#'            attribute of a [get_fs()] result passed as `data` (which may
#'            be absent — `merMod` results carry none); an explicit `fsb`
#'            always wins.
#' @param reduce Controls how per-unit `fsL`/`fsT`/`fsb` from a group fitted
#'            with missing data (`missing = "fiml"`; per-pattern lists of
#'            matrices), per-cluster values from a `merMod` fit (3-D arrays),
#'            or per-observation values from a `mirt` fit, are collapsed to a
#'            single representative value per group for
#'            stage 2. A no-op when the per-unit quantities are constant within
#'            the group (e.g. complete single-group data). A group (or the
#'            whole data) with no scorable rows at all (every row missing) is
#'            an error: there is nothing to pool. One of `"mean"`
#'            (the default) or `"median"`. With `"mean"` the pooled `fsT` is a
#'            convex combination of the per-unit (positive semi-definite)
#'            matrices and so remains positive semi-definite; with `"median"`
#'            the reduction is element-wise and need not be, in which case a
#'            warning is emitted when the pooled `fsT` is not positive
#'            semi-definite. The pooled (not the nested/per-cluster) `fsT` and
#'            `fsL` are what get attached to the returned fit.
#' @param vfsLT The sampling covariance matrix of the free `fsL`/`fsT`
#'            elements, taken from the `vfsLT` attribute of a [get_fs()]
#'            result fitted with `vfsLT = TRUE`. Required when
#'            `corrected_se = TRUE`; ignored otherwise.
#' @param corrected_se A logical; when `TRUE`, the stage-2 covariance of the
#'            returned fit is replaced by the first-order (delta-method)
#'            correction of [vcov_corrected()] and the `tspa_corrected`
#'            attribute is set to `TRUE`. Requires a multi-factor fit (both
#'            `fsT` and `fsL` supplied) and `vfsLT`. Supported for
#'            single-group and multigroup fits (the multigroup `fsL`/`fsT`
#'            are the per-group list attributes from `get_fs()`). Default
#'            `FALSE` (the returned fit is unchanged).
#' @param which_free An optional numeric vector of positions selecting which
#'            `fsL`/`fsT` free elements to propagate through the corrected
#'            covariance (see [vcov_corrected()]); used only when
#'            `corrected_se = TRUE`.
#' @param product A logical; when `TRUE`, [tspa()] automatically computes
#'            the double-mean-centered product indicators (via
#'            [compute_fs_prod()]) for every latent in `model` that names
#'            the product of two of the model's factor scores — by
#'            concatenation (the latent `xm` is the product of the scores
#'            of `x` and `m`) or in lavaan's interaction syntax (`x:m`,
#'            rendered under the concatenated name) — and incorporates
#'            them into the stage-2 measurement model, so
#'            `tspa("y ~ x + m + x:m", data = get_fs(...), product = TRUE)`
#'            needs no pre-computed product columns. In the single-factor
#'            (score-scale) path the product latent loads 1 on its
#'            indicator with error variance the (per-group pooled,
#'            `reduce`) product SE, like every other single-factor latent;
#'            in the multi-factor path it loads with the implied loading
#'            `gamma` and error variance `se_P^2`, both evaluated at the
#'            (pooled) `fsL`/`fsT` with the `psi` attribute. When two
#'            product latents share a factor score (e.g. `x:m` and `x:z`
#'            share `x`), their indicators' measurement errors are
#'            correlated, and `tspa()` fixes those error covariances in the
#'            stage-2 model, estimated from the stage-1 `fsL`/`fsT`/`psi`
#'            (the Isserlis expansion of the joint-normal score-error
#'            moments — the same joint-normality assumptions the
#'            [compute_fs_prod()] product SE rests on); with a single
#'            product latent there is nothing to fix. Single-group
#'            models only (v1); not supported with `corrected_se = TRUE`.
#'            Default `FALSE`, which leaves the manual workflow
#'            ([get_fs()] with `product` set, or [compute_fs_prod()], up
#'            front, the product SE in `se_fs`) unchanged.
#' @param ... Additional arguments passed to \code{\link[lavaan]{sem}}. See
#'            \code{\link[lavaan]{lavOptions}} for a complete list.
#' @return An object of class \code{lavaan} carrying the following
#'         \code{R2spa}-specific attributes: \code{tspaModel}, the stage-2
#'         model syntax actually fitted; \code{tspa_call}, the matched
#'         `tspa()` call; and \code{tspa_args}, the argument list used to
#'         build the stage-2 model (captured at fit time as evaluated
#'         values), which lets the fit be re-evaluated without its original
#'         environment, e.g. by the `vcov_corrected()` SE correction. On a
#'         fit whose measurement inputs were derived from `data` (see
#'         `Details`), \code{tspa_args} carries the *resolved* values
#'         (post-pooling plain matrices, the derived `se_fs`), so a
#'         replay via `do.call(tspa, attr(fit, "tspa_args"))` re-passes
#'         them explicitly, skips the derivation, and cannot double-pool.
#'         When
#'         \code{fsT}/\code{fsL} are supplied (multi-factor measurement
#'         model), the (possibly reduced) matrices are also attached as the
#'         \code{fsT}/\code{fsL} attributes, and \code{pooled_fs} records the
#'         \code{reduce} method used when per-unit values were collapsed.
#'         When \code{corrected_se = TRUE}, the returned fit additionally
#'         carries \code{tspa_corrected = TRUE} and its covariance is the
#'         first-order corrected matrix, so `vcov()`, `se()`, and
#'         `standardizedSolution()` on it report the corrected standard
#'         errors, for multigroup fits as well as single-group ones.
#'
#' @export
#'
#' @examples
#' library(lavaan)
#'
#' # single-group, two-factor example, factor scores obtained separately
#' # get factor scores
#' fs_dat_ind60 <- get_fs(object = PoliticalDemocracy,
#'                        model = "ind60 =~ x1 + x2 + x3")
#' fs_dat_dem60 <- get_fs(object = PoliticalDemocracy,
#'                        model = "dem60 =~ y1 + y2 + y3 + y4")
#' fs_dat <- cbind(fs_dat_ind60, fs_dat_dem60)
#' # tspa model
#' tspa(model = "dem60 ~ ind60", data = fs_dat,
#'      se_fs = c(ind60 = fs_dat_ind60[1, "fs_ind60_se"],
#'                dem60 = fs_dat_dem60[1, "fs_dem60_se"]))
#'
#' # the same fit with se_fs derived from the fs_<v>_se columns of the
#' # cbind()ed get_fs() results (no explicit se_fs)
#' tspa(model = "dem60 ~ ind60", data = fs_dat)
#'
#' # single-group, three-factor example
#' mod2 <- "
#'   # latent variables
#'     ind60 =~ x1 + x2 + x3
#'     dem60 =~ y1 + y2 + y3 + y4
#'     dem65 =~ y5 + y6 + y7 + y8
#' "
#' fs_dat2 <- get_fs(PoliticalDemocracy, model = mod2, std.lv = TRUE)
#' tspa(model = "dem60 ~ ind60
#'               dem65 ~ ind60 + dem60",
#'      data = fs_dat2,
#'      fsT = attr(fs_dat2, "fsT"),
#'      fsL = attr(fs_dat2, "fsL"))
#'
#' # the same fit with the measurement inputs derived from the get_fs()
#' # result (no explicit fsT/fsL)
#' tspa(model = "dem60 ~ ind60
#'               dem65 ~ ind60 + dem60",
#'      data = fs_dat2)
#'
#' # multigroup, two-factor example
#' mod3 <- "
#'   # latent variables
#'     visual =~ x1 + x2 + x3
#'     speed =~ x7 + x8 + x9
#' "
#' fs_dat3 <- get_fs(HolzingerSwineford1939, model = mod3, std.lv = TRUE,
#'                   group = "school")
#' tspa(model = "visual ~ speed",
#'      data = fs_dat3,
#'      fsT = attr(fs_dat3, "fsT"),
#'      fsL = attr(fs_dat3, "fsL"),
#'      group = "school")
#'
#' # the same fit with the measurement inputs derived from the get_fs()
#' # result (no explicit fsT/fsL)
#' tspa(model = "visual ~ speed",
#'      data = fs_dat3,
#'      group = "school")
#'
#' # multigroup, three-factor example
#' mod4 <- "
#'   # latent variables
#'     visual =~ x1 + x2 + x3
#'     textual =~ x4 + x5 + x6
#'     speed =~ x7 + x8 + x9
#' "
#' fs_dat4 <- get_fs(HolzingerSwineford1939, model = mod4, std.lv = TRUE,
#'                   group = "school")
#' tspa(model = "visual ~ speed
#'               textual ~ visual + speed",
#'      data = fs_dat4,
#'      fsT = attr(fs_dat4, "fsT"),
#'      fsL = attr(fs_dat4, "fsL"),
#'      group = "school")
#'
#' # get factor scores
#' fs_dat_visual <- get_fs(object = HolzingerSwineford1939,
#'                         model = "visual =~ x1 + x2 + x3",
#'                         group = "school",
#'                         format = "list")
#' fs_dat_speed <- get_fs(object = HolzingerSwineford1939,
#'                        model = "speed =~ x7 + x8 + x9",
#'                        group = "school",
#'                        format = "list")
#' fs_hs <- cbind(do.call(rbind, fs_dat_visual),
#'                do.call(rbind, fs_dat_speed))
#'
#' # tspa model
#' tspa(model = "visual ~ speed",
#'      data = fs_hs,
#'      se_fs = data.frame(visual = c(0.3391326, 0.311828),
#'                         speed = c(0.2786875, 0.2740507)),
#'      group = "school",
#'      group.equal = "regressions")
#'
#' # manually adding equality constraints on the regression coefficients
#' tspa(model = "visual ~ c(b1, b1) * speed",
#'      data = fs_hs,
#'      se_fs = list(visual = c(0.3391326, 0.311828),
#'                   speed = c(0.2786875, 0.2740507)),
#'      group = "school")
#'
#' # Missing data (FIML): per-pattern fsL/fsT are pooled within the group
#' data("HolzingerSwineford1939", package = "lavaan")
#' hs <- HolzingerSwineford1939
#' set.seed(1334)
#' hs$x2[!rbinom(nrow(hs), 1, 0.4)] <- NA
#' hs$x8[!rbinom(nrow(hs), 1, 0.4)] <- NA
#' mod_fin <- "
#'   visual =~ x1 + x2 + x3
#'   speed  =~ x7 + x8 + x9
#' "
#' fit_fin <- suppressWarnings(cfa(mod_fin, data = hs, missing = "fiml"))
#' fs_fin <- get_fs(fit_fin)
#' tspa("visual ~ speed", data = fs_fin,
#'      fsT = attr(fs_fin, "fsT"), fsL = attr(fs_fin, "fsL"),
#'      reduce = "mean")
#' # opt-in element-wise reduction (may lose positive semi-definiteness)
#' suppressWarnings(tspa("visual ~ speed", data = fs_fin,
#'      fsT = attr(fs_fin, "fsT"), fsL = attr(fs_fin, "fsL"),
#'      reduce = "median"))
#'
#' # merMod: per-cluster fsL/fsT are pooled (one value per cluster)
#' library(lme4)
#' lmod <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' fs_mer <- get_fs(lmod)
#' tspa("u1 ~ u0", data = fs_mer,
#'      fsT = attr(fs_mer, "fsT"), fsL = attr(fs_mer, "fsL"))
#'
#' # Product-score auto-compute (opt-in): the model's `x:m` latent is the
#' # product of the `x` and `m` scores (rendered under the concatenated
#' # name `xm`), computed on the fly from the data's stage-1 attributes
#' # (no pre-computed product columns needed)
#' set.seed(2116)
#' covx <- matrix(c(1, 0.4, 0.4, 1), 2)
#' eta <- as.data.frame(MASS::mvrnorm(500, rep(0, 2), covx))
#' names(eta) <- c("x", "m")
#' lk <- list(x = c(0.9, 0.8, 0.7), m = c(0.85, 0.75, 0.65),
#'            y = c(0.75, 0.7, 0.65))
#' etay <- 0.5 * eta$x + 0.4 * eta$m + 0.3 * eta$x * eta$m
#' obs <- setNames(lapply(c("x", "m"), function(v0) {
#'   eta[[v0]] %*% t(lk[[v0]]) + rnorm(1500)
#' }), c("x", "m"))
#' obs$y <- etay %*% t(lk$y) + rnorm(1500)
#' df <- as.data.frame(do.call(cbind, obs))
#' names(df) <- c(paste0("x", 1:3), paste0("m", 1:3), paste0("y", 1:3))
#' fs_prod <- get_fs(df, model = "x =~ x1 + x2 + x3
#'                     m =~ m1 + m2 + m3
#'                     y =~ y1 + y2 + y3", std.lv = TRUE)
#' tspa("y ~ x + m + x:m", data = fs_prod, product = TRUE)


tspa <- function(model, data, reliability = NULL, se = "standard",
                 se_fs = NULL, fsT = NULL, fsL = NULL, fsb = NULL,
                 reduce = c("mean", "median"),
                 vfsLT = NULL, corrected_se = FALSE, which_free = NULL,
                 product = FALSE,
                 ...) {
  reduce <- match.arg(reduce)
  # Set when the per-unit fsT/fsL/fsb attributes were collapsed to a single
  # representative per-group set via `reduce` (PLAN 09); attached to the
  # returned fit as the "pooled_fs" attribute.
  pooled_fs <- NULL

  if (!inherits(model, "character")) {
    stop("The structural path model provided is not a string.")
  }

  if (!is.null(reliability)) {
    stop("tspa() currently does not support reliability model")
  }
  if (!is.character(se)) {
    warning("using `se` to set se for factor scores is deprecated. ",
            "use `se_fs` instead.")
  }
  if (!is.logical(product) || length(product) != 1L || is.na(product)) {
    stop("'product' must be a single TRUE/FALSE value.", call. = FALSE)
  }
  if (isTRUE(corrected_se) && isTRUE(product)) {
    stop(
      "'corrected_se = TRUE' is not supported with 'product = TRUE' ",
      "(v1): the product indicator's stage-1 uncertainty is not ",
      "propagated by the delta-method correction.",
      call. = FALSE
    )
  }

  # PLAN 13: capture before the coercion below (a missing `se_fs` becomes a
  # 0 x 0 data frame there, so NULL-ness is the "the user supplied se_fs"
  # signal the measurement-input derivations gate on).
  se_fs_given <- !is.null(se_fs)

  if (!is.data.frame(se_fs)) {
    se_fs <- as.data.frame(as.list(se_fs))
  }
  if (ncol(se_fs) > 0L &&
      any(vapply(se_fs, function(x) !is.numeric(x), logical(1)))) {
    stop("'se_fs' must contain numeric standard errors.", call. = FALSE)
  }
  if (xor(is.null(fsT), is.null(fsL))) {
    stop("Please provide both or none of fsT and fsL.")
  }

  # PLAN 13: derive the measurement inputs from a get_fs() result when the
  # caller omitted them. Explicit arguments always win (D2); a supplied
  # se_fs suppresses the multi-factor derivation (D3); with both forms
  # available and nothing passed, the multi-factor (attribute) form wins
  # (D4). Derivation fires only for argument values that are NULL, and every
  # such call errors today (no measurement inputs reach stage 2), so no
  # currently-working call changes behavior.
  derived_prov_err <- NULL
  if (is.null(fsT) && is.null(fsL) && !se_fs_given) {
    # Multi-factor derivation: the data's own fsT/fsL (and fsb) attributes
    # (unified or list format; per-unit shapes -- FIML per-pattern, merMod
    # per-cluster, mirt per-obs -- are reduced by the pooling below).
    attr_T <- attr(data, "fsT")
    attr_L <- attr(data, "fsL")
    if (!is.null(attr_T) && !is.null(attr_L)) {
      # Provenance gate: derive only from attributes that resolve as a
      # get_fs() result, reusing resolve_fs_per_row()'s informative errors
      # instead of duplicating its shape rules. A hand-rolled frame with
      # plain matrix attributes but no get_fs() provenance is NOT derived;
      # it falls through to the fail-fast error below, which carries the
      # gate's message.
      prov <- tryCatch(
        {
          resolve_fs_per_row(data)
          TRUE
        },
        error = function(e) conditionMessage(e)
      )
      if (isTRUE(prov)) {
        fsT <- attr_T
        fsL <- attr_L
        fsb <- attr(data, "fsb") # may be NULL (merMod has none)
      } else {
        derived_prov_err <- prov
      }
    }
    # Single-factor derivation: the data's own fs_<v>/fs_<v>_se columns
    # (a cbind'd get_fs() result carries no attributes, so its se_fs is
    # derived from the columns).
    if (is.null(fsT)) {
      derived_se <- derive_sf_se_fs(data, list(...)[["group"]], reduce)
      if (!is.null(derived_se)) {
        se_fs <- derived_se
      }
    }
  }
  # Fail fast with an actionable message instead of lavaan's "model is NULL
  # or not a valid type for it!" when no measurement inputs exist.
  if (is.null(fsT) && nrow(se_fs) == 0) {
    stop(
      "No measurement inputs found for the factor scores in 'data'. ",
      "Please supply one of: (1) 'se_fs' (single-factor), (2) 'fsT' and ",
      "'fsL' (multi-factor), or (3) a get_fs() result as 'data' (its ",
      "attributes, or its fs_<v>/fs_<v>_se columns, carry the inputs).",
      if (is.null(derived_prov_err)) {
        ""
      } else {
        paste0(
          " The 'fsT'/'fsL' attributes on 'data' were detected but not ",
          "used because the data does not look like a get_fs() result: ",
          derived_prov_err
        )
      },
      call. = FALSE
    )
  }

  # A plain matrix stands for a single group, so a length-1 list may be mixed
  # with a plain matrix (single-group unified get_fs() output); only differing
  # group counts are an error.
  if (!is.null(fsT)) {
    nT <- if (is.list(fsT)) length(fsT) else 1
    nL <- if (is.list(fsL)) length(fsL) else 1
    if (nT != nL) {
      stop(
        if (nT > nL) {
          "'fsL' must be a list of the same length as 'fsT' for a multigroup model."
        } else {
          "'fsT' must be a list of the same length as 'fsL' for a multigroup model."
        }
      )
    }
    # Per-unit values -- FIML per-pattern (per-group attribute lists of
    # matrices), merMod per-cluster (3-D arrays), or mirt per-obs (the
    # `mirt_per_obs` marker; PLAN 11) -- are collapsed to a single
    # representative per-group fsT/fsL/fsb via `reduce` (PLAN 09/11) BEFORE
    # multigroup detection, name matching, and schema building, so those
    # run on clean per-group (or single) matrices.
    if (is_per_unit_fs(fsT, fsL,
                       mirt_per_obs = isTRUE(attr(data, "mirt_per_obs")) ||
                         isTRUE(attr(data, "per_obs")))) {
      pooled <- pool_per_unit(data, reduce, have_int = !is.null(fsb))
      fsT <- pooled$fsT
      fsL <- pooled$fsL
      if (!is.null(fsb)) {
        fsb <- pooled$fsb
      }
      pooled_fs <- reduce
    } else if (!is.null(fsb) && any(vapply(
      # Residual backstop: a per-pattern fsb (list of vectors) without
      # matching per-unit fsT/fsL is still unsupported and, unwrapped, would
      # be misread as a multigroup vector list.
      if (is.list(fsb)) fsb else list(fsb),
      function(e) is.list(e) && length(e) > 0L &&
        all(
          vapply(e, function(v) is.vector(v) && !is.matrix(v), logical(1))
        ),
      logical(1)
    ))) {
      stop(
        "tspa() does not yet support groups with multiple missing-data ",
        "patterns: the 'fsb' attribute contains a per-pattern list of ",
        "vectors without matching per-unit 'fsT'/'fsL'. Fit the stage-1 ",
        "model on complete data.",
        call. = FALSE
      )
    }
  }
  multigroup <- if (!is.null(fsT)) {
    is.list(fsT) && length(fsT) > 1
  } else {
    nrow(se_fs) > 1
  }

  if (!is.null(fsT)) {
    if (multigroup) {
      fs_names <- colnames(fsT[[1]])
      dat_names <- if (is.data.frame(data)) names(data) else names(data[[1]])
    } else {
      fs_names <- if (is.list(fsT)) colnames(fsT[[1]]) else colnames(fsT)
      dat_names <- names(data)
    }
    names_match <- fs_names %in% dat_names
    if (any(!names_match)) {
      stop(
        "Names of factor score variables do not match those in the input data."
      )
    }
    # The schema reads the score names from fsT's column names and the
    # indicator names from fsL's row names; missing dimnames would render
    # NA names into the model string (a cryptic rbind failure downstream,
    # not an error here).
    L1n <- if (is.list(fsL)) fsL[[1L]] else fsL
    T1n <- if (is.list(fsT)) fsT[[1L]] else fsT
    if (is.null(rownames(T1n)) || is.null(colnames(T1n)) ||
        is.null(rownames(L1n)) || is.null(colnames(L1n))) {
      stop(
        "'fsT' and 'fsL' must carry both row and column names (the ",
        "factor-score names); the stage-2 model is built from them.",
        call. = FALSE
      )
    }
  }

  if (multigroup && is.null(list(...)[["group"]])) {
    stop("Please specify 'group = ' to fit a multigroup model in lavaan.")
  }

  if (is.null(fsT)) { # single-factor measurement model
    # FIML per-pattern se pooling (PLAN 09): a get_fs() result fitted with
    # missing data carries per-observation `fs_<v>_se` columns that vary
    # within a group, which a single (per-group) `se_fs` row cannot
    # represent. Detect within-group se variation and, when it fires,
    # replace `se_fs` with the per-group reduction of the per-row se
    # columns. A group whose se column is constant within the group gives
    # no signal, so complete-data (homogeneous) behavior is unchanged, and
    # data without the se columns (not a get_fs() result) is left as-is.
    # A list of per-group data frames (a get_fs(format = "list") result)
    # is coerced to a single frame first: lavaan's data= must be a
    # data.frame. The coercion lives here rather than at the top of
    # tspa() because rbind drops the data's custom attributes (group_col,
    # and the per-group fsT/fsL/fsb the multi-factor per-unit pooling
    # above needs un-rbind'd); on this branch fsT is NULL, so per-unit
    # pooling cannot apply. The group column is resolved with the
    # derivation's convention (attribute, else the group= argument, else
    # a literal "group" column); the attribute is read pre-rbind, since
    # rbind drops it.
    sf_group_col <- attr(data, "group_col")
    if (inherits(data, "list")) {
      data <- do.call(rbind, data)
    }
    if (is.null(sf_group_col)) {
      g_arg <- list(...)[["group"]]
      if (is.character(g_arg) && length(g_arg) == 1L &&
          g_arg %in% names(data)) {
        sf_group_col <- g_arg
      } else if ("group" %in% names(data)) {
        sf_group_col <- "group"
      }
    }
    if (!is.null(sf_group_col) && sf_group_col %in% names(data) &&
        !is.null(se_fs) && nrow(se_fs) > 0 && ncol(se_fs) > 0) {
      se_cols <- paste0("fs_", colnames(se_fs), "_se")
      if (all(se_cols %in% names(data))) {
        varied <- any(vapply(se_cols, function(cl) {
          any(vapply(
            split(data[[cl]], data[[sf_group_col]]),
            function(x) sum(!is.na(unique(x))) > 1L,
            logical(1)
          ))
        }, logical(1)))
        if (varied) {
          se_fs <- pool_se_fs(data, colnames(se_fs), reduce, sf_group_col)
        }
      }
    }
    # Pooling may have grown se_fs to one row per group: refresh the
    # multigroup flag (computed before this branch) and the group=
    # requirement, which must track the value actually fitted.
    multigroup <- nrow(se_fs) > 1
    if (multigroup && is.null(list(...)[["group"]])) {
      stop("Please specify 'group = ' to fit a multigroup model in lavaan.")
    }
    # Opt-in product-score auto-compute: for each model latent that names
    # the product of two of this model's factor scores — concatenated
    # (`xm` for `x` and `m`) or in lavaan interaction syntax (`x:m`,
    # rendered under the concatenated name) — ensure the DMC product
    # columns exist and the product SE enters se_fs. The score-scale
    # convention of this path is unchanged: the product latent loads 1 on
    # its indicator, error = the (pooled) product SE, like every other
    # single-factor latent.
    prods <- NULL
    if (isTRUE(product)) {
      prods <- tspa_product_latents(model, names(data), colnames(se_fs))
      if (!is.null(prods)) {
        if (multigroup) {
          stop(
            "'product = TRUE' is not supported for multigroup models ",
            "(v1: single-group only).",
            call. = FALSE
          )
        }
        model <- tspa_rewrite_product_toks(model, prods)
        # An explicit product SE keyed by the model token (`x:m`) is
        # renamed to the render name (`xm`) the generated model uses.
        # The token column may arrive as the check.names() form (`x.m`,
        # from the as.data.frame() coercion or data.frame()'s default) or
        # as the literal token (a data.frame built with check.names =
        # FALSE).
        for (k in seq_len(nrow(prods))) {
          if (prods$tok[k] == prods$v[k]) next
          tk <- make.names(prods$tok[k])
          for (nm in unique(c(tk, prods$tok[k]))) {
            if (nm %in% colnames(se_fs)) {
              se_fs[[prods$v[k]]] <- se_fs[[nm]]
              se_fs[[nm]] <- NULL
              break
            }
          }
        }
        data <- tspa_ensure_product_cols(data, prods)
        for (k in seq_len(nrow(prods))) {
          v <- prods$v[k]
          if (v %in% colnames(se_fs)) next
          se_col <- tspa_product_se_col(data, prods$a[k], prods$b[k])
          se_v <- if (nrow(se_fs) > 1L) {
            pool_se_col(data, se_col, reduce, sf_group_col)
          } else {
            pool_se_col(data, se_col, reduce)
          }
          se_fs[[v]] <- se_v
        }
      }
    }
    # Product-score columns (compute_fs_prod: `fs_a:fs_b`) are not valid
    # lavaan variable names; the schema's generated model name for latent
    # `v` is `fs_v`, so a matching product-score column is aliased into a
    # working copy of the data. Manual pre-renames keep working (the alias
    # is a no-op when `fs_v` already exists).
    al <- tspa_sf_alias(data, se_fs)
    data <- al$data
    # Product-indicator error covariances: with more than one product
    # latent, two products sharing a factor score have correlated
    # measurement errors that the stage-2 model must fix. The authoritative
    # pair set is `prods` (product = TRUE); the alias's `prod_map` covers
    # the manual workflow (pre-computed product columns with the product
    # SEs in se_fs), including replayed/pre-renamed fits where the alias
    # itself was a no-op.
    prod_pairs <- prods
    if (is.null(prod_pairs) && !is.null(al$prod_map)) {
      prod_pairs <- al$prod_map
    }
    prod_ecov <- NULL
    if (!is.null(prod_pairs) && nrow(prod_pairs) >= 2L) {
      fsL_a <- attr(data, "fsL")
      fsT_a <- attr(data, "fsT")
      psi_a <- attr(data, "psi")
      if (is.null(fsL_a) || is.null(fsT_a) || is.null(psi_a)) {
        stop(
          "Cannot compute the product-indicator error covariances: the ",
          "data lacks the stage-1 attributes (fsL/fsT/psi) that a direct ",
          "single-group get_fs() result carries (a cbind()ed result drops ",
          "them). Pass the un-cbind()ed get_fs() result.",
          call. = FALSE
        )
      }
      L1 <- if (is.list(fsL_a)) fsL_a[[1L]] else fsL_a
      T1 <- if (is.list(fsT_a)) fsT_a[[1L]] else fsT_a
      psi1 <- fs_psi_matrix(psi_a)
      prod_ecov <- tspa_prod_ecov(prod_pairs, L1, T1, psi1)
      # The schema emits the ecov rows for group 1 only (single-group
      # v1); a multigroup fit would leave the other groups' product
      # indicators' correlated errors unmodeled. A get_fs() multigroup
      # result already errors above (its per-group psi list is rejected
      # by fs_psi_matrix()); this catches hand-built plain-matrix
      # attributes.
      if (multigroup && !is.null(prod_ecov)) {
        stop(
          "The product-indicator error covariances are not supported for ",
          "multigroup models (v1: single-group only).",
          call. = FALSE
        )
      }
    }
    # A non-finite SE (an all-NA fs_<v>_se group, a degenerate product SE
    # computation, or user input) would become a fixed NaN in the stage-2
    # model: fail with an actionable message instead of lavaan's parse
    # failure on the NaN. (Column-wise: is.finite() has no data.frame
    # method.)
    if (any(vapply(se_fs, function(cl) !all(is.finite(cl)), logical(1)))) {
      stop(
        "The 'se_fs' values fed to the stage-2 model are not all finite ",
        "(check the 'fs_<v>_se' columns of 'data' and any product SE ",
        "computation).",
        call. = FALSE
      )
    }
    tspaModel <- tspa_sf(model, data, se_fs, prod_ecov)
  } else { # multi-factor measurement model
    # Opt-in product-score auto-compute, multi-factor path: the product
    # latent (named by concatenation or lavaan interaction syntax, as in
    # the single-factor path) loads on its DMC indicator with the implied
    # loading `gamma` (the true-latent scale of this path) and error
    # variance `se_P^2`, both evaluated at the (pooled) fsL/fsT with the
    # psi attribute.
    prods_mf <- NULL
    prods_ecov <- NULL
    if (isTRUE(product)) {
      if (is.list(fsT) && length(fsT) > 1L) {
        stop(
          "'product = TRUE' is not supported for multigroup models ",
          "(v1: single-group only).",
          call. = FALSE
        )
      }
      L1 <- if (is.list(fsL)) fsL[[1L]] else fsL
      dat_names <- if (is.data.frame(data)) names(data) else names(data[[1L]])
      prods <- tspa_product_latents(model, dat_names, colnames(L1))
      if (!is.null(prods)) {
        if (!is.data.frame(data)) {
          stop(
            "'product = TRUE' requires a single-group data-frame 'data'; ",
            "a list of per-group data frames is not supported (v1).",
            call. = FALSE
          )
        }
        model <- tspa_rewrite_product_toks(model, prods)
        T1 <- if (is.list(fsT)) fsT[[1L]] else fsT
        psi <- fs_psi_matrix(attr(data, "psi"))
        data <- tspa_ensure_product_cols(data, prods)
        ld_vals <- numeric(nrow(prods))
        se2_vals <- numeric(nrow(prods))
        for (k in seq_len(nrow(prods))) {
          i <- match(prods$a[k], colnames(L1))
          j <- match(prods$b[k], colnames(L1))
          ld_vals[k] <- fs_prod_gamma(L1, i, j)
          se2_vals[k] <- fs_prod_se2(L1, T1, psi, i, j)
        }
        # A non-positive se_P^2 is a degenerate stage-1 (the single-factor
        # path's sqrt_or_na() maps it to NA and the finiteness guard there
        # errors); a fixed negative error variance would otherwise be
        # silently fitted.
        bad_se2 <- !is.finite(se2_vals) | se2_vals < 0
        if (any(bad_se2)) {
          stop(
            "The implied product error variance (se_P^2) for the product ",
            "latent(s) ", paste(prods$v[bad_se2], collapse = ", "),
            " is not positive and finite; check the stage-1 fit.",
            call. = FALSE
          )
        }
        prods_mf <- data.frame(
          v = prods$v, ld = ld_vals, se2 = se2_vals,
          stringsAsFactors = FALSE
        )
        for (k in seq_len(nrow(prods))) {
          tgt <- paste0("fs_", prods$v[k])
          src <- tspa_product_col(data, prods$a[k], prods$b[k])
          if (is.na(src)) {
            stop(
              "The model references the product latent '", prods$v[k],
              "' but the data has no matching product-score column.",
              call. = FALSE
            )
          }
          if (!(tgt %in% names(data))) {
            data[[tgt]] <- data[[src]]
          }
        }
        # Product-indicator error covariances: pairs of product latents
        # sharing a factor score have correlated measurement errors that
        # the stage-2 model fixes (NULL for a single product, as for a
        # pair set with no nonzero covariances).
        prods_ecov <- tspa_prod_ecov(prods, L1, T1, psi)
      }
    }
    # A non-finite fsT/fsL/fsb (hand-built input; the pooled inputs are
    # already guarded) would be rendered as "NA" in the model string,
    # which lavaan parses as a parameter label: the fixed measurement
    # value would silently become a free estimate.
    if (!all_finite_values(fsT) || !all_finite_values(fsL) ||
        (!is.null(fsb) && !all_finite_values(fsb))) {
      stop(
        "'fsT', 'fsL', or 'fsb' contains non-finite values; the stage-2 ",
        "model cannot be built from them (an 'NA' fixed value would be ",
        "parsed by lavaan as a parameter label, silently unfixing it).",
        call. = FALSE
      )
    }
    tspaModel <- tspa_mf(model, data, fsT, fsL, fsb, prods_mf, prods_ecov)
    if (inherits(data, "list")) {
      data <- do.call(rbind, data)
    }
  }

  tspa_fit <- sem(model = tspaModel,
                  data  = data,
                  se = se,
                  ...)
  attr(tspa_fit, "tspaModel") <- tspaModel
  if (!is.null(fsT)) {
    attr(tspa_fit, "fsT") <- fsT
    attr(tspa_fit, "fsL") <- fsL
    if (!is.null(pooled_fs)) {
      attr(tspa_fit, "pooled_fs") <- pooled_fs
    }
  }
  # Self-contained replay record: every refit-relevant argument as an
  # evaluated value, plus the forwarded lavaan dots, so
  # do.call(tspa, attr(fit, "tspa_args")) reproduces the fit without the
  # original call environment. All element names are tspa() formals, so the
  # spliced dots structurally cannot shadow them.
  attr(tspa_fit, "tspa_args") <- c(
    list(model = model, data = data, reliability = reliability, se = se,
         se_fs = se_fs, fsT = fsT, fsL = fsL, fsb = fsb, reduce = reduce,
         vfsLT = vfsLT, corrected_se = corrected_se,
         which_free = which_free, product = product),
    list(...)
  )
  attr(tspa_fit, "tspa_call") <- match.call()
  # First-order (delta-method) SE correction: replace the lavaan covariance
  # with vcov(fit) + J %*% vfsLT %*% t(J) (see vcov_corrected()). Only
  # meaningful for the multi-factor path (fsT and fsL supplied) with a
  # vfsLT matrix from get_fs(..., vfsLT = TRUE); single- and multi-group
  # fits.
  if (isTRUE(corrected_se)) {
    if (is.null(fsT) || is.null(fsL) || is.null(vfsLT)) {
      stop(
        "corrected_se = TRUE requires a multi-factor fit (both 'fsT' and ",
        "'fsL' supplied) and a 'vfsLT' matrix from get_fs(..., vfsLT = TRUE)."
      )
    }
    corrected <- vcov_corrected(tspa_fit, vfsLT = vfsLT,
                                which_free = which_free)
    # The @vcov slot write goes through the lavaan-internal boundary
    # (single point of coupling, canary-covered).
    tspa_fit <- tsp_set_vcov(tspa_fit, corrected)
    attr(tspa_fit, "tspa_corrected") <- TRUE
  }
  return(tspa_fit)
}

# ---------------------------------------------------------------------------
# Per-unit factor-score pooling (PLAN 09): collapse FIML per-pattern or
# merMod per-cluster fsL/fsT/fsb attributes to a single representative
# per-group set that the existing stage-2 machinery accepts. The
# long-form expansion reuses resolve_fs_per_row()/fs_row_cols() (R/fs_indiv.R)
# so the per-row values are exactly the ones fs_indiv() reports.
# ---------------------------------------------------------------------------

# TRUE when `fsT`/`fsL` carry per-unit heterogeneity and are poolable:
# a 3-D array (merMod per-cluster), a per-group value that is itself a
# list of matrices (lavaan per-pattern under missing data), or mirt
# per-obs input (a flat list of bare matrices, one per observation) when
# the authoritative `mirt_per_obs` marker is passed through. The marker is
# required so a bare flat list of matrices is never mistaken for per-unit
# values on its own (PLAN 09 Section 8 / PLAN 11).
is_per_unit_fs <- function(fsT, fsL, mirt_per_obs = FALSE) {
  per_unit <- function(x) {
    (is.array(x) && length(dim(x)) == 3L) ||
      (is.list(x) && any(vapply(
        x,
        function(e) is.list(e) && length(e) > 0L &&
          all(vapply(e, is.matrix, logical(1))),
        logical(1)
      )))
  }
  per_unit(fsT) || per_unit(fsL) ||
    (mirt_per_obs && is.list(fsT) && length(fsT) > 1L)
}

# Pool one get_fs() result (`fs`) to a single representative fsT/fsL/fsb
# per group. `reduce` is "mean" (default; the mean of PSD matrices is PSD)
# or "median" (element-wise; may break PSD, guarded by a warning). Returns
# list(fsT, fsL, fsb): for a multi-group result (group_vals non-null) the
# reduction is done within each group and each component is a list of one
# value per group, named by group label in the group column's
# first-appearance order (the lavaan stage-2 group order, verified for
# character and factor columns; for lavaan/merMod input this equals the
# stage-1 attribute list order, for mirt input the stage-1 list is in
# mirt's level order, which can differ); otherwise each component is a
# single matrix/vector (SG FIML, merMod). The returned shapes are exactly
# what the existing stage-2 schema accepts.
pool_per_unit <- function(fs, reduce, have_int) {
  resolved <- resolve_fs_per_row(fs)
  # Effective per-row grouping. lavaan/merMod carry it in the resolved
  # structure; mirt per-obs does not (resolve_per_obs() always sets the
  # group fields to NULL), so recover it from the data's own `group` column
  # when the mirt marker is set -- otherwise MG mirt would pool into one
  # group (PLAN 11).
  g_vals <- resolved$group_vals
  g_col <- resolved$group_col
  mirt_mg <- (isTRUE(attr(fs, "mirt_per_obs")) ||
               isTRUE(attr(fs, "per_obs"))) &&
    is.data.frame(fs) && "group" %in% names(fs)
  if (is.null(g_vals) && mirt_mg) {
    g_col <- "group"
    g_vals <- as.character(fs[["group"]])
  }
  ref_T <- resolved$blocks[[1L]]$fsT
  ref_L <- resolved$blocks[[1L]]$fsL
  q <- ncol(ref_T)
  has_int <- have_int && !is.null(resolved$blocks[[1L]]$fsb)

  n <- resolved$n
  k_ld <- q * q
  k_ev <- q * (q + 1L) / 2
  k_int <- if (has_int) q else 0L

  se_mat <- matrix(NA_real_, nrow = n, ncol = q)
  ld_mat <- matrix(NA_real_, nrow = n, ncol = k_ld)
  ev_mat <- matrix(NA_real_, nrow = n, ncol = k_ev)
  int_mat <- if (has_int) matrix(NA_real_, nrow = n, ncol = k_int) else NULL

  # Same block loop as fs_indiv(): expand every unit's per-pattern/
  # per-cluster matrices once per member row (the per-observation-equal
  # reading), na.rm in the reduction below drops the all-NA rows. A block
  # no row maps to (inconsistent hand-built input) is skipped rather than
  # warning on the empty assignment.
  for (b in seq_along(resolved$blocks)) {
    blk <- resolved$blocks[[b]]
    rows_b <- which(resolved$pattern_idx == b)
    if (length(rows_b) == 0L) next
    vals <- fs_row_cols(
      resolved$scores[rows_b, , drop = FALSE],
      blk$fsL,
      blk$fsT,
      if (has_int) blk$fsb else NULL
    )
    se_mat[rows_b, ] <- vals[, seq_len(q), drop = FALSE]
    ld_mat[rows_b, ] <- vals[, (q + 1L):(q + k_ld), drop = FALSE]
    ev_mat[rows_b, ] <- vals[, (q + k_ld + 1L):(q + k_ld + k_ev), drop = FALSE]
    if (has_int) {
      int_mat[rows_b, ] <-
        vals[, (q + k_ld + k_ev + 1L):ncol(vals), drop = FALSE]
    }
  }

  reduce_fn <- if (reduce == "median") stats::median else mean
  dn_L <- dimnames(ref_L)
  dn_T <- dimnames(ref_T)

  # Reduce one row subset to (fsL, fsT, fsb): ld column-major into a
  # q x q matrix; ev row-major lower triangle (i-outer / j<=i-inner, the
  # fs_row_cols order) back into a symmetric matrix; int named from the
  # score rows.
  pool_rows <- function(rows_i) {
    ld_red <- vapply(seq_len(k_ld), function(k) {
      reduce_fn(ld_mat[rows_i, k], na.rm = TRUE)
    }, numeric(1))
    ev_red <- vapply(seq_len(k_ev), function(k) {
      reduce_fn(ev_mat[rows_i, k], na.rm = TRUE)
    }, numeric(1))
    int_red <- if (has_int) {
      vapply(seq_len(k_int), function(k) {
        reduce_fn(int_mat[rows_i, k], na.rm = TRUE)
      }, numeric(1))
    } else {
      NULL
    }
    fsL_p <- matrix(ld_red, nrow = q, ncol = q)
    dimnames(fsL_p) <- dn_L
    fsT_p <- matrix(NA_real_, nrow = q, ncol = q)
    count <- 1L
    for (i in seq_len(q)) {
      for (j in seq_len(i)) {
        fsT_p[i, j] <- ev_red[count]
        fsT_p[j, i] <- ev_red[count]
        count <- count + 1L
      }
    }
    dimnames(fsT_p) <- dn_T
    fsb_p <- if (is.null(int_red)) NULL else {
      names(int_red) <- rownames(ref_T)
      int_red
    }
    list(fsT = fsT_p, fsL = fsL_p, fsb = fsb_p)
  }
  psd_guard <- function(res, label) {
    # A reduction over no scorable rows (a group whose rows are all
    # missing, or an unused level of the group column) yields an all-NA
    # pooled fsT: fail with an actionable message instead of the
    # misleading PSD warning (which would suggest switching reduce while
    # the real problem is the absence of data).
    if (all(!is.finite(res$fsT))) {
      stop(
        "tspa() pooling found no scorable rows in ",
        if (is.null(label)) "the data" else paste0("group '", label, "'"),
        " (all rows are missing or unscorable), so the pooled 'fsT' is ",
        "empty. Check the missing data and the group column.",
        call. = FALSE
      )
    }
    emin <- pooled_fsT_min_eigen(res$fsT)
    if (!is.finite(emin) || emin < -.Machine$double.eps^0.5) {
      warning(
        "Pooled 'fsT'",
        if (is.null(label)) "" else paste0(" for group '", label, "'"),
        " (reduce = \"", reduce, "\") is not positive semi-definite; ",
        "consider reduce = \"mean\" (the mean of PSD matrices is PSD).",
        call. = FALSE
      )
    }
  }

  if (!is.null(g_vals)) {
    # Multi-group: reduce within each group. Each component (fsT/fsL/fsb)
    # is a list of one value per group, named by group label in the data's
    # group order (factor: level order, the lavaan stage-2 group order),
    # matching the stage-1 attribute list order.
    glabs <- if (!is.null(g_col) &&
                is.data.frame(fs) &&
                g_col %in% names(fs)) {
      # First-appearance order of the data's own group column: verified
      # on lavaan 0.7-2 to be lavaan's `group =` stage-2 order for BOTH
      # character and factor columns (a factor's LEVEL order is not used
      # by lavaan). For a factor this also drops NA (the completely-
      # missing rows of a mirt result), which pool_rows() would
      # otherwise turn into an empty "group"; every mirt group has at
      # least one scorable row, so all K groups appear exactly once.
      # (A mirt stage-1 attribute list is in mirt's groupNames/level
      # order; the pooled list must match the stage-2 order instead,
      # which differs when the data rows are not in level order.)
      gc <- fs[[g_col]]
      if (is.factor(gc)) unique(gc[!is.na(gc)]) else unique(gc)
    } else {
      # list-format input (no group column in the named list): group order
      # of the list, which the per-row values follow.
      unique(g_vals)
    }
    T_list <- vector("list", length(glabs))
    L_list <- vector("list", length(glabs))
    b_list <- vector("list", length(glabs))
    gnames <- as.character(glabs)
    names(T_list) <- gnames
    names(L_list) <- gnames
    names(b_list) <- gnames
    for (k in seq_along(glabs)) {
      rows_g <- which(g_vals == gnames[k])
      res <- pool_rows(rows_g)
      psd_guard(res, gnames[k])
      T_list[[k]] <- res$fsT
      L_list[[k]] <- res$fsL
      b_list[[k]] <- res$fsb
    }
    return(list(
      fsT = T_list,
      fsL = L_list,
      fsb = if (has_int) b_list else NULL
    ))
  }
  res <- pool_rows(seq_len(n))
  psd_guard(res, NULL)
  res
}

# Smallest eigenvalue of a symmetric pooled fsT (non-finite when the
# matrix is non-finite, e.g. a group with no scorable rows).
pooled_fsT_min_eigen <- function(T_mat) {
  if (!all(is.finite(T_mat))) {
    return(NA_real_)
  }
  min(eigen(T_mat, symmetric = TRUE, only.values = TRUE)$values)
}

# TRUE when every value of a matrix/vector (or of each element of a list of
# them, the per-group form) is finite. Non-finite fixed values would be
# rendered as "NA"/"NaN" in the model string, which lavaan parses as a
# parameter label: the fixed value would silently become a free estimate.
all_finite_values <- function(x) {
  all(vapply(if (is.list(x)) x else list(x),
             function(e) !is.null(e) && all(is.finite(e)), logical(1)))
}

# Per-group reduction of one materialized `<se_col>` column of a get_fs()
# result: a scalar (no group column), or one value per group in the
# stage-2 group order. `rows` restricts the reduction to a row subset
# (pool_se_fs() pools each group through it).
pool_se_col <- function(data, se_col, reduce, group_col = NULL, rows = NULL) {
  reduce_fn <- if (reduce == "median") stats::median else mean
  x <- data[[se_col]]
  if (!is.null(rows)) {
    return(reduce_fn(x[rows], na.rm = TRUE))
  }
  if (is.null(group_col)) {
    return(reduce_fn(x, na.rm = TRUE))
  }
  gvals <- data[[group_col]]
  gnames <- fs_group_order(gvals)
  # unname: vapply names the result with a character input's values
  unname(vapply(gnames, function(g) reduce_fn(x[gvals == g], na.rm = TRUE),
                numeric(1)))
}

# Per-group reduction of the materialized `fs_<v>_se` columns of a get_fs()
# result (single-factor FIML): one row per group (group order of the
# data's group column) with the per-latent pooled se; a numeric vector
# without a group column.
pool_se_fs <- function(data, se_names, reduce, group_col) {
  if (is.null(group_col)) {
    return(vapply(se_names, function(v) {
      pool_se_col(data, paste0("fs_", v, "_se"), reduce)
    }, numeric(1)))
  }
  gvals <- data[[group_col]]
  gnames <- fs_group_order(gvals)
  out <- do.call(rbind, lapply(gnames, function(g) {
    vals <- vapply(se_names, function(v) {
      pool_se_col(data, paste0("fs_", v, "_se"), reduce,
                  rows = (gvals == g))
    }, numeric(1))
    data.frame(t(vals), check.names = FALSE)
  }))
  colnames(out) <- se_names
  rownames(out) <- gnames
  out
}

# ---------------------------------------------------------------------------
# Product-score auto-compute (opt-in `product = TRUE`): the product latents
# a model string names, the missing product columns, and the per-pair
# column lookups shared by the single- and multi-factor paths.
# ---------------------------------------------------------------------------

# The product latents a model string names. An identifier token of the
# model (a data column or a numeric value is never one) is a product latent
# in either form:
#   (a) lavaan interaction syntax `a:b` — a single token with non-empty
#       parts, both distinct known factor scores; it renders under the
#       concatenated name `ab`, because in the generated model `a:b` would
#       be parsed by lavaan as an interaction of the (latent) variables;
#   (b) the concatenation of two distinct known factor scores (e.g. `xm`
#       for `x` and `m`).
# Tokens that are neither (a label, an unknown name, an `a:b` token whose
# parts are not both factor scores — passed through to lavaan as an
# ordinary interaction — ...) are ignored: a genuinely unknown model
# variable still fails later in lavaan, as before product support. A known
# factor score that is also a concatenation is a product latent only when
# the data has no score column for it (a present `fs_v` column means the
# regular latent wins and there is nothing to compute). A candidate
# matching two DIFFERENT pairs is ambiguous and errors, as is a model that
# names the same pair twice (`x:m` and `xm`) or whose render name collides
# with another model variable. Returns a data frame (tok, v, a, b) — `tok`
# the model token as written, `v` the render name — or NULL when no
# product latent is named.
tspa_product_latents <- function(model, data_names, known) {
  mtxt <- sub("(?m)#.*$", "", paste(model, collapse = "\n"), perl = TRUE)
  # `:` is kept inside tokens (the interaction syntax `a:b`); everything
  # else non-name splits.
  toks <- unique(unlist(strsplit(mtxt, "[^A-Za-z0-9_.:]+")))
  toks <- toks[nzchar(toks)]
  cand <- setdiff(toks, data_names)
  cand <- cand[!grepl("^[+-]?[0-9]+(\\.[0-9]*)?([eE][+-]?[0-9]+)?$", cand)]
  known <- unique(known)
  prods <- list()
  for (v in cand) {
    # Form (a): `a:b`.
    m <- regmatches(
      v, regexec("^([A-Za-z0-9_.]+):([A-Za-z0-9_.]+)$", v, perl = TRUE)
    )[[1L]]
    if (length(m) == 3L) {
      a <- m[2L]
      b <- m[3L]
      if (a != b && a %in% known && b %in% known) {
        prods[[length(prods) + 1L]] <-
          c(tok = v, v = paste0(a, b), a = a, b = b)
      }
      next
    }
    # Form (b): the concatenation of two known factor scores. A known
    # factor score with this name wins over the product reading whenever
    # the data carries its score column.
    if (v %in% known && paste0("fs_", v) %in% data_names) next
    hits <- list()
    for (i in seq_along(known)) {
      for (j in seq_along(known)) {
        if (i == j) next
        if (paste0(known[i], known[j]) == v) {
          hits[[length(hits) + 1L]] <- c(known[i], known[j])
        }
      }
    }
    if (length(hits) == 0L) next
    keys <- vapply(hits, function(p) paste(p, collapse = ":"), character(1))
    if (length(unique(keys)) > 1L) {
      stop(
        "Cannot determine which factor-score pair the model variable '",
        v, "' is the product of (",
        paste(unique(keys), collapse = ", "),
        "). Rename it to disambiguate.",
        call. = FALSE
      )
    }
    p1 <- hits[[1L]]
    prods[[length(prods) + 1L]] <- c(tok = v, v = v, a = p1[1L], b = p1[2L])
  }
  if (length(prods) == 0L) {
    return(NULL)
  }
  df <- data.frame(
    tok = vapply(prods, function(p) p["tok"], character(1)),
    v = vapply(prods, function(p) p["v"], character(1)),
    a = vapply(prods, function(p) p["a"], character(1)),
    b = vapply(prods, function(p) p["b"], character(1)),
    stringsAsFactors = FALSE
  )
  # The same pair named twice (`x:m` and `xm`) would render two model
  # variables with the same name.
  pair_keys <- vapply(
    seq_len(nrow(df)),
    function(k) paste(sort(c(df$a[k], df$b[k])), collapse = ":"),
    character(1)
  )
  if (any(duplicated(pair_keys))) {
    dups <- df$tok[duplicated(pair_keys) |
                      duplicated(pair_keys, fromLast = TRUE)]
    stop(
      "The model names the same factor-score pair more than once (",
      paste(unique(dups), collapse = ", "),
      "). Name it once.",
      call. = FALSE
    )
  }
  # A render name colliding with another model variable (e.g. a stage-1
  # latent called `xm` while the model also names `x:m`) would silently
  # merge two latents.
  for (k in seq_len(nrow(df))) {
    if (df$v[k] != df$tok[k] &&
        df$v[k] %in% setdiff(toks, c(df$tok[k], data_names))) {
      stop(
        "The render name '", df$v[k], "' of the product latent '",
        df$tok[k], "' collides with another variable of the model. ",
        "Rename the pair.",
        call. = FALSE
      )
    }
  }
  df
}

# Rewrite the model string's product tokens to their render names (form
# (a) tokens only; form (b) tokens already are their render names). The
# token is matched as a whole model variable: its boundaries are
# non-name characters.
tspa_rewrite_product_toks <- function(model, prods) {
  for (k in seq_len(nrow(prods))) {
    if (identical(prods$tok[k], prods$v[k])) next
    esc <- gsub(".", "\\\\.", prods$tok[k], fixed = TRUE)
    model <- gsub(
      paste0("(?<![A-Za-z0-9_.])", esc, "(?![A-Za-z0-9_.])"),
      prods$v[k], model, perl = TRUE
    )
  }
  model
}

# The existing product-score column of the pair (a, b) in the data, in
# EITHER orientation (`fs_a:fs_b` or `fs_b:fs_a`; a user may have computed
# either), NA when absent.
tspa_product_col <- function(data, a, b) {
  cands <- c(paste0("fs_", a, ":fs_", b), paste0("fs_", b, ":fs_", a))
  hit <- intersect(cands, names(data))
  if (length(hit) == 0L) NA_character_ else hit[1L]
}

# The existing product-SE column of the pair (a, b), NA when absent.
tspa_product_se_col <- function(data, a, b) {
  cl <- tspa_product_col(data, a, b)
  if (is.na(cl)) return(NA_character_)
  paste0(cl, "_se")
}

# Ensure the DMC product columns for the pairs named by `prods` exist in
# `data`: a pair missing its column in both orientations is computed via
# compute_fs_prod() (which validates the input and rejects multigroup,
# merMod, per-observation, and attribute-less data informatively). The
# new columns are created in canonical (sorted) pair order. Data without
# the stage-1 attributes (e.g. a cbind()ed get_fs() result, which drops
# them) cannot be auto-computed and is rejected with the remedy.
tspa_ensure_product_cols <- function(data, prods) {
  pairs <- unique(prods[, c("a", "b"), drop = FALSE])
  missing <- vapply(
    seq_len(nrow(pairs)),
    function(k) is.na(tspa_product_col(data, pairs$a[k], pairs$b[k])),
    logical(1)
  )
  if (any(missing)) {
    if (is.null(attr(data, "fsL")) || is.null(attr(data, "fsT")) ||
        is.null(attr(data, "psi"))) {
      stop(
        "Cannot auto-compute the product column(s) for ",
        paste(vapply(
          seq_len(nrow(pairs))[missing],
          function(k) paste(sort(c(pairs$a[k], pairs$b[k])),
                            collapse = ":"),
          character(1)
        ), collapse = ", "),
        ": the data lacks the stage-1 attributes (fsL/fsT/psi) that a ",
        "direct single-group get_fs() result carries (a cbind()ed result ",
        "drops them). Pre-compute the product columns with get_fs(product = ) ",
        "or compute_fs_prod(), or pass the un-cbind()ed get_fs() result.",
        call. = FALSE
      )
    }
    spec <- paste(
      vapply(
        seq_len(nrow(pairs))[missing],
        function(k) paste(sort(c(pairs$a[k], pairs$b[k])), collapse = ":"),
        character(1)
      ),
      collapse = " + "
    )
    data <- compute_fs_prod(data, product = spec)
  }
  data
}

# Product-indicator error covariances for the stage-2 model: for every
# unordered pair of the product latents named by `prods` (a data frame with
# columns v, a, b — `a`/`b` the factor scores of the DMC indicator
# `fs_a:fs_b`; the extra `tok` column of tspa_product_latents() is ignored),
# the measurement-error covariance of the two indicators,
# fs_prod_ecov() (R/compute_fs_prod.R), evaluated at the (pooled) matrices
# L/T/psi. Two product indicators sharing a factor score (e.g. `x:m` and
# `x:z` share `x`) have correlated measurement errors that the stage-2 model
# must fix; a (numerically) zero value — the two pairs linked by no
# score-error moment — is dropped. Returns a data frame (v1, v2, ecov) in
# (r, p) loop order, or NULL when there are fewer than two products or every
# pair is zero.
tspa_prod_ecov <- function(prods, L, Tm, psi) {
  if (nrow(prods) < 2L) {
    return(NULL)
  }
  if (!is.matrix(L) || !is.matrix(Tm) || !is.matrix(psi)) {
    stop(
      "The product-indicator error covariances require plain-matrix ",
      "stage-1 'fsL'/'fsT'/'psi' values; per-pattern (FIML) blocks are not ",
      "supported with more than one product latent (v1).",
      call. = FALSE
    )
  }
  ia <- match(prods$a, colnames(L))
  ib <- match(prods$b, colnames(L))
  if (anyNA(ia) || anyNA(ib)) {
    stop(
      "The product-indicator error covariances cannot be computed: the ",
      "factor score(s) ",
      paste(sort(unique(c(prods$a[is.na(ia)], prods$b[is.na(ib)]))),
            collapse = ", "),
      " are not a column of the stage-1 'fsL' matrix.",
      call. = FALSE
    )
  }
  n <- nrow(prods)
  v1 <- character()
  v2 <- character()
  ecov <- numeric()
  for (r in seq_len(n - 1L)) {
    for (p in seq(r + 1L, n)) {
      val <- fs_prod_ecov(L, Tm, psi, ia[p], ib[p], ia[r], ib[r])
      if (abs(val) > 1e-12) {
        v1 <- c(v1, prods$v[r])
        v2 <- c(v2, prods$v[p])
        ecov <- c(ecov, val)
      }
    }
  }
  if (length(ecov) == 0L) {
    return(NULL)
  }
  data.frame(v1 = v1, v2 = v2, ecov = ecov, stringsAsFactors = FALSE)
}

# Stage-2 group order (PLAN 13): the unique values of the group column in
# first-appearance order -- verified on lavaan 0.7-2 to be lavaan's own
# `group =` order for BOTH character and factor columns (lavaan does not
# use a factor's level order; `unique()` on a factor also yields
# first-appearance, not level, order). The same convention the stage-1
# attribute list order and pool_per_unit() follow. Shared by pool_se_fs()
# and the derived se_fs, which must be row-aligned with the stage-2 groups.
fs_group_order <- function(gvals) {
  as.character(unique(gvals))
}

# Single-factor derivation (PLAN 13): build the `se_fs` an attribute-less
# get_fs() result (e.g. a cbind'd one) implies from its own columns. The
# latent set is every simple score column `fs_<v>` that carries a matching
# numeric `fs_<v>_se` column (product-score columns `fs_a:fs_b` are out of
# scope, D7), in first-appearance order of the score columns (the cbind
# order of the canonical vignette flows). With no group column a single
# row is returned; with one, one row per group in the stage-2 group order
# with the within-group se reduced by `reduce` (a no-op for the constant
# columns of complete data, idempotent with the FIML se pooling that may
# run on the result downstream). The group column is the data's `group_col`
# attribute, else the user's `group` argument (a cbind'd frame has no
# attribute, so the argument is the only stage-2 group signal), else a
# literal "group" column. Returns NULL when no latent pair exists (the
# caller keeps the empty se_fs).
derive_sf_se_fs <- function(data, group_arg, reduce) {
  dnames <- names(data)
  if (is.null(dnames)) {
    return(NULL)
  }
  score_cols <- dnames[startsWith(dnames, "fs_") & !grepl(":", dnames)]
  latents <- character()
  for (cl in score_cols) {
    v <- substr(cl, 4L, nchar(cl))
    se_cl <- paste0("fs_", v, "_se")
    if (nzchar(v) && se_cl %in% dnames && is.numeric(data[[se_cl]])) {
      latents <- c(latents, v)
    }
  }
  if (length(latents) == 0L) {
    return(NULL)
  }
  grp_col <- attr(data, "group_col")
  if (!is.null(grp_col) && !(grp_col %in% dnames)) {
    grp_col <- NULL
  }
  if (is.null(grp_col)) {
    if (is.character(group_arg) && length(group_arg) == 1L &&
        group_arg %in% dnames) {
      grp_col <- group_arg
    } else if ("group" %in% dnames) {
      grp_col <- "group"
    }
  }
  pooled <- pool_se_fs(data, latents, reduce, grp_col)
  if (is.data.frame(pooled)) {
    return(pooled)
  }
  # No group column: pool_se_fs() returned a named numeric vector; shape
  # it the way the se_fs coercion does (one row, one column per latent).
  as.data.frame(as.list(pooled))
}

# ---------------------------------------------------------------------------
# Stage-2 model schema (PLAN 04): owned and frozen by R2spa, not by lavaan.
# One row per (statement term, group):
#   lhs   statement left-hand side (NA for verbatim user rows)
#   op    "=~" | "~~" | "~" | "raw" ("raw" rows carry user syntax)
#   rhs   right-hand variable ("1" for intercepts); "raw" rows carry the
#         user's line verbatim
#   value per-group fixed value (NA for "raw" rows)
#   free  0 = fixed (every injected row); NA for "raw" rows
#   group 1-based group index (NA for "raw" rows)
#   label R2spa-generated label for the statement (__r2spa_ldN__ /
#         __r2spa_evN__ / __r2spa_intN__): a generated namespace user
#         labels/variable names cannot collide with; carried as row
#         metadata, never emitted into the syntax
#   kind  "user" | "struct" | "error_var" | "error_cov" | "intercept"
# ---------------------------------------------------------------------------

tspa_row <- function(lhs, op, rhs, value, group, kind, label) {
  data.frame(
    lhs = if (is.null(lhs)) NA_character_ else lhs,
    op = op,
    rhs = rhs,
    value = value,
    free = if (op == "raw") NA_integer_ else 0L,
    group = if (is.null(group)) NA_integer_ else as.integer(group),
    label = label,
    kind = kind,
    stringsAsFactors = FALSE
  )
}

# The user model string, carried verbatim as a single "raw" row (no
# re-parsing and no round-trip through line splitting, so leading/trailing
# newlines survive exactly as written).
tspa_user_rows <- function(model) {
  tspa_row(NA, "raw", paste0(model, collapse = "\n"), NA, NA, "user", NA)
}

# Single-factor (se_fs) schema: per latent, one fixed-loading struct row
# and one error-variance row per group; values follow the se_fs rows
# (groups) in order; fixed error-covariance rows between product
# indicators sharing a factor score when `prod_ecov` is given
# (single-group v1, hence group 1).
tspa_schema_sf <- function(model, se, prod_ecov = NULL) {
  var <- colnames(se)
  fs <- paste0("fs_", var)
  ng <- nrow(se)
  rows <- list(tspa_user_rows(model))
  for (k in seq_along(var)) {
    lab <- paste0("__r2spa_ld", k, "__")
    for (g in seq_len(ng)) {
      rows[[length(rows) + 1L]] <- tspa_row(
        var[k], "=~", fs[k], 1, g, "struct", lab
      )
    }
    ev_lab <- paste0("__r2spa_ev", k, "__")
    for (g in seq_len(ng)) {
      rows[[length(rows) + 1L]] <- tspa_row(
        fs[k], "~~", fs[k], se[g, k]^2, g, "error_var", ev_lab
      )
    }
  }
  if (!is.null(prod_ecov)) {
    for (k in seq_len(nrow(prod_ecov))) {
      rows[[length(rows) + 1L]] <- tspa_row(
        paste0("fs_", prod_ecov$v1[k]), "~~", paste0("fs_", prod_ecov$v2[k]),
        prod_ecov$ecov[k], 1L, "error_cov",
        paste0("__r2spa_peck", k, "__")
      )
    }
  }
  do.call(rbind, rows)
}

# Multi-factor (fsT/fsL/fsb) schema: per latent, one struct row per score
# term and per group; error rows follow the lower triangle (incl. diagonal)
# of fsT in column-major order — the legacy per-group value routing made
# explicit and unit-testable; per-score intercept rows when fsb is given;
# product-indicator rows when `prods` is given (one fixed loading row and
# one fixed error-variance row per product latent, single-group v1) and
# fixed error-covariance rows between product indicators sharing a factor
# score when `prod_ecov` is given (group 1 likewise).
tspa_schema_mf <- function(model, fsT, fsL, fsb, prods = NULL,
                           prod_ecov = NULL) {
  # `fsT`/`fsL` are plain matrices for a single-group model or named lists
  # of them for a multigroup model. Single-group unified get_fs() output
  # carries length-1 list attributes, so either shape is accepted on either
  # side (e.g. list-valued `fsT` with a plain identity `fsL` for Bartlett).
  if (is.list(fsT) && length(fsT) > 1) {
    if (!is.list(fsL) || length(fsL) != length(fsT)) {
      stop("'fsL' must be a list of the same length as 'fsT' for a ",
           "multigroup model.")
    }
    ngroup <- length(fsT)
    L_list <- fsL
    T_list <- fsT
  } else if (is.list(fsL) && length(fsL) > 1) {
    if (!is.list(fsT) || length(fsT) != length(fsL)) {
      stop("'fsT' must be a list of the same length as 'fsL' for a ",
           "multigroup model.")
    }
    ngroup <- length(fsL)
    L_list <- fsL
    T_list <- fsT
  } else {
    ngroup <- 1
    L_list <- if (is.list(fsL)) fsL else list(fsL)
    T_list <- if (is.list(fsT)) fsT else list(fsT)
  }
  fsL1 <- L_list[[1]]
  var <- colnames(fsL1)
  fs <- rownames(fsL1)

  rows <- list(tspa_user_rows(model))
  for (k in seq_along(var)) {
    lab <- paste0("__r2spa_ld", k, "__")
    for (i in seq_along(fs)) {
      for (g in seq_len(ngroup)) {
        rows[[length(rows) + 1L]] <- tspa_row(
          var[k], "=~", fs[i], L_list[[g]][i, k], g, "struct", lab
        )
      }
    }
  }
  ev_count <- 0L
  tri <- which(!upper.tri(T_list[[1]]), arr.ind = TRUE)
  for (k in seq_len(nrow(tri))) {
    i <- tri[k, 1]
    j <- tri[k, 2]
    ev_count <- ev_count + 1L
    lab <- paste0("__r2spa_ev", ev_count, "__")
    kind <- if (i == j) "error_var" else "error_cov"
    for (g in seq_len(ngroup)) {
      rows[[length(rows) + 1L]] <- tspa_row(
        fs[i], "~~", fs[j], T_list[[g]][i, j], g, kind, lab
      )
    }
  }
  if (!is.null(fsb)) {
    B_list <- if (is.list(fsb)) fsb else list(fsb)
    for (i in seq_along(fs)) {
      lab <- paste0("__r2spa_int", i, "__")
      for (g in seq_len(ngroup)) {
        rows[[length(rows) + 1L]] <- tspa_row(
          fs[i], "~", "1", B_list[[g]][i], g, "intercept", lab
        )
      }
    }
  }
  # Product indicators: the DMC product column `fs_v` (aliased from
  # `fs_a:fs_b` upstream) is a single fixed indicator of the product
  # latent — loading `gamma`, error variance `se_P^2`, evaluated at the
  # (pooled) matrices upstream (single-group v1, hence group 1); two
  # product indicators sharing a factor score get a fixed error-covariance
  # row (from `prod_ecov`, likewise group 1).
  if (!is.null(prods)) {
    for (k in seq_len(nrow(prods))) {
      fv <- paste0("fs_", prods$v[k])
      rows[[length(rows) + 1L]] <- tspa_row(
        prods$v[k], "=~", fv, prods$ld[k], 1L, "struct",
        paste0("__r2spa_pld", k, "__")
      )
      rows[[length(rows) + 1L]] <- tspa_row(
        fv, "~~", fv, prods$se2[k], 1L, "error_var",
        paste0("__r2spa_pev", k, "__")
      )
    }
    if (!is.null(prod_ecov)) {
      for (k in seq_len(nrow(prod_ecov))) {
        rows[[length(rows) + 1L]] <- tspa_row(
          paste0("fs_", prod_ecov$v1[k]), "~~",
          paste0("fs_", prod_ecov$v2[k]),
          prod_ecov$ecov[k], 1L, "error_cov",
          paste0("__r2spa_peck", k, "__")
        )
      }
    }
  }
  do.call(rbind, rows)
}

# Statements of one schema: rows with the given identity grouped into
# consecutive statements (order of first appearance).
tspa_statements <- function(sch, id) {
  id <- match(id, unique(id))
  starts <- c(1L, which(diff(id) != 0L) + 1L)
  lapply(starts, function(i) sch[id == id[i], , drop = FALSE])
}

# Per-group values of a single-term statement, in group order.
tspa_stmt_values <- function(st) {
  gs <- sort(unique(st$group))
  unlist(lapply(gs, function(g) st$value[st$group == g]))
}

# Value string for one statement term: a bare number for a single (i.e.
# single-group) value, `c(v1, v2, ...)` for a multigroup term. Empirically
# verified on lavaan 0.7.2: the bare forms (loading, fixed self-variance,
# intercept, and a negative value, e.g. `a + -0.02 * b`) parse identically
# to the `c(...)` forms -- same row order, estimates, and vcov() -- while a
# bare minus between terms (`a - 0.02 * b`) is rejected by lavaan. Terms are
# therefore always joined with " + " and the sign is carried by the value.
tspa_vals_str <- function(vals) {
  if (length(vals) == 1L) {
    as.character(vals[[1L]])
  } else {
    paste0("c(", paste0(vals, collapse = ", "), ")")
  }
}

# Per-term value string for a possibly multi-term statement; `terms` are the
# ordered unique rhs values (one row value per group each).
tspa_stmt_vals <- function(st, terms) {
  trm <- match(st$rhs, terms)
  paste(
    vapply(seq_along(terms), function(k) {
      tspa_vals_str(tspa_stmt_values(st[trm == k, , drop = FALSE]))
    }, character(1)),
    collapse = " + "
  )
}

# The single renderer (PLAN 04): schema -> lavaan model syntax string.
# Every statement is emitted as `lhs <sp> op <sp> rhs` (spaces around `*`
# kept) with the per-term values bare for single-group statements and
# `c(v1, v2, ...)` for multigroup ones; the user model block is carried
# verbatim. The emitted string is semantically identical to the legacy
# c()-everywhere builder (row order, estimates, and vcov() unchanged).
tspa_render <- function(sch, style = c("sf", "mf")) {
  style <- match.arg(style)
  user_lines <- sch$rhs[sch$kind == "user"]
  struct <- sch[sch$kind == "struct", , drop = FALSE]
  errors <- sch[sch$kind %in% c("error_var", "error_cov"), , drop = FALSE]
  ints <- sch[sch$kind == "intercept", , drop = FALSE]
  if (style == "sf") {
    latent_var_str <- paste(
      vapply(tspa_statements(struct, struct$lhs), function(st) {
        paste0(st$lhs[1], " =~ ", tspa_stmt_vals(st, st$rhs[1]),
               " * ", st$rhs[1], "\n")
      }, character(1)),
      collapse = ""
    )
    error_constraint_str <- paste(
      vapply(tspa_statements(errors, paste(errors$lhs, errors$rhs,
                                           sep = "|")),
             function(st) {
               paste0(st$lhs[1], " ~~ ", tspa_stmt_vals(st, st$rhs[1]),
                      " * ", st$rhs[1], "\n")
             }, character(1)),
      collapse = ""
    )
    elems <- c(
      "# latent variables (indicated by factor scores)",
      latent_var_str,
      "# constrain the errors",
      error_constraint_str,
      "# structural model",
      paste(user_lines, collapse = "\n")
    )
  } else {
    latent_var_str <- vapply(
      tspa_statements(struct, struct$lhs),
      function(st) {
        terms <- unique(st$rhs)
        loadings_c <- paste(
          vapply(terms, function(t) {
            tr <- st[st$rhs == t, , drop = FALSE]
            paste0(tspa_stmt_vals(tr, t), " * ", t)
          }, character(1)),
          collapse = " + "
        )
        paste0("# latent variables (indicated by factor scores)\n",
               st$lhs[1], " =~ ", loadings_c)
      },
      character(1)
    )
    error_constraint_str <- vapply(
      tspa_statements(errors, paste(errors$lhs, errors$rhs, sep = "|")),
      function(st) {
        paste0("# constrain the errors\n", st$lhs[1], " ~~ ",
               tspa_stmt_vals(st, st$rhs[1]), " * ", st$rhs[1])
      },
      character(1)
    )
    if (nrow(ints) > 0) {
      intercept_constraint <- vapply(
        tspa_statements(ints, ints$lhs),
        function(st) {
          paste0("# constrain the intercepts\n", st$lhs[1], " ~ ",
                 tspa_stmt_vals(st, st$rhs[1]), " * 1")
        },
        character(1)
      )
    } else {
      # The legacy builder emitted an empty element here (rendered as a
      # blank line before the structural block) when no intercepts exist.
      intercept_constraint <- ""
    }
    elems <- c(
      latent_var_str,
      error_constraint_str,
      intercept_constraint,
      "# structural model",
      paste(user_lines, collapse = "\n")
    )
  }
  paste0(elems, collapse = "\n")
}

tspa_sf <- function(model, data, se = NULL, prod_ecov = NULL) {
  if (!is.null(se) && nrow(se) > 0L) {
    return(tspa_render(tspa_schema_sf(model, se, prod_ecov), style = "sf"))
  }
}

tspa_mf <- function(model, data, fsT, fsL, fsb, prods = NULL,
                    prod_ecov = NULL) {
  tspa_render(tspa_schema_mf(model, fsT, fsL, fsb, prods, prod_ecov),
              style = "mf")
}

# ---------------------------------------------------------------------------
# Product-score (compute_fs_prod) auto-alias: the schema's generated model
# name for latent `v` is `fs_v`; a data column `fs_a:fs_b` (a,b latent
# names in se_fs) whose names concatenate to `v` is copied into the working
# data as `fs_v`. Old user models that pre-rename the product-score column
# keep working because the alias is a no-op when `fs_v` already exists.
# The matched (a, b) pair of every se entry `v` with a single matching
# product-score column is also returned as `prod_map` (columns v, a, b;
# NULL when none) — including the pre-renamed/replayed cases where the
# alias itself is a no-op — because the product-indicator error
# covariances (tspa_prod_ecov()) need the pair even then.
# ---------------------------------------------------------------------------

tspa_sf_alias <- function(data, se) {
  is_lst <- inherits(data, "list") && !is.data.frame(data)
  dnames <- if (is_lst) names(data[[1]]) else names(data)
  se_names <- colnames(se)
  aliases <- character()
  pm_v <- character()
  pm_a <- character()
  pm_b <- character()
  for (v in se_names) {
    tgt <- paste0("fs_", v)
    cand <- character()
    cand_ab <- character()
    for (col in dnames) {
      pos <- regexpr(":", col, fixed = TRUE)
      if (pos < 1L) next
      a <- substr(col, 1L, pos - 1L)
      b <- substr(col, pos + 1L, nchar(col))
      # Product-score columns are `fs_a:fs_b` (both parts score names).
      if (!grepl("^fs_", a) || !grepl("^fs_", b)) next
      a <- sub("^fs_", "", a)
      b <- sub("^fs_", "", b)
      if (!(a %in% se_names) || !(b %in% se_names)) next
      if (paste0(a, b) == v || paste0(b, a) == v) {
        cand <- c(cand, col)
        cand_ab <- c(cand_ab, a, b)
      }
    }
    if (length(cand) == 0) next
    if (length(cand) > 1) {
      if (tgt %in% dnames) next # pre-renamed: the no-op keeps working
      stop(
        "Cannot determine which product-score column in the input data ",
        "corresponds to the latent variable '", v, "': ",
        paste0("\"", cand, "\"", collapse = " or "),
        ". Rename it to \"", tgt, "\" to disambiguate."
      )
    }
    if (tgt %in% dnames) {
      # The score column already exists (a user pre-rename, or a replay of
      # a previously aliased fit): no alias, but `v` is still the product
      # of (a, b) — the error covariances need the pair.
    } else {
      if (is_lst) {
        for (i in seq_along(data)) data[[i]][[tgt]] <- data[[i]][[cand]]
      } else {
        data[[tgt]] <- data[[cand]]
      }
      aliases <- c(aliases, paste0(tgt, " <- ", cand))
    }
    pm_v <- c(pm_v, v)
    pm_a <- c(pm_a, cand_ab[1L])
    pm_b <- c(pm_b, cand_ab[2L])
  }
  prod_map <- if (length(pm_v) > 0L) {
    data.frame(v = pm_v, a = pm_a, b = pm_b, stringsAsFactors = FALSE)
  } else {
    NULL
  }
  list(data = data, aliases = aliases, prod_map = prod_map)
}
