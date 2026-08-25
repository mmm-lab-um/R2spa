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
#'            the group (e.g. complete single-group data). One of `"mean"`
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


tspa <- function(model, data, reliability = NULL, se = "standard",
                 se_fs = NULL, fsT = NULL, fsL = NULL, fsb = NULL,
                 reduce = c("mean", "median"),
                 vfsLT = NULL, corrected_se = FALSE, which_free = NULL,
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

  # PLAN 13: capture before the coercion below (a missing `se_fs` becomes a
  # 0 x 0 data frame there, so NULL-ness is the "the user supplied se_fs"
  # signal the measurement-input derivations gate on).
  se_fs_given <- !is.null(se_fs)

  if (!is.data.frame(se_fs)) {
    se_fs <- as.data.frame(as.list(se_fs))
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
    names_match <- lapply(fs_names, function(x) x %in% dat_names) |> unlist()
    if (any(!names_match)) {
      stop(
        "Names of factor score variables do not match those in the input data."
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
    sf_group_col <- attr(data, "group_col")
    if (is.null(sf_group_col) && "group" %in% names(data)) {
      sf_group_col <- "group"
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
    # Product-score columns (get_fs_int: `fs_a:fs_b`) are not valid lavaan
    # variable names; the schema's generated model name for latent `v` is
    # `fs_v`, so a matching product-score column is aliased into a working
    # copy of the data. Manual pre-renames keep working (the alias is a
    # no-op when `fs_v` already exists).
    data <- tspa_sf_alias(data, se_fs)$data
    tspaModel <- tspa_sf(model, data, se_fs)
  } else { # multi-factor measurement model
    tspaModel <- tspa_mf(model, data, fsT, fsL, fsb)
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
         which_free = which_free),
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
# value per group, named by group label (in the data's group order -- the
# same order the stage-1 attribute lists use); otherwise each component is
# a single matrix/vector (SG FIML, merMod). The returned shapes are exactly
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
  # reading), na.rm in the reduction below drops the all-NA rows.
  for (b in seq_along(resolved$blocks)) {
    blk <- resolved$blocks[[b]]
    rows_b <- which(resolved$pattern_idx == b)
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
      # Unique values of the data's own group column: a factor's level
      # order (the lavaan stage-2 group order), a character vector's
      # first-appearance order -- matching the stage-1 attribute list
      # order. mirt's `group` is a factor; use the full level set (not
      # unique()) so the pooled list always has K entries, one per mirt
      # group, matching lavaan's `group=` levels.
      gc <- fs[[g_col]]
      if (isTRUE(mirt_mg) && is.factor(gc)) levels(gc) else unique(gc)
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

# Per-group reduction of the materialized `fs_<v>_se` columns of a get_fs()
# result (single-factor FIML): one row per group (group order of the
# data's group column) with the per-latent pooled se; a numeric vector
# without a group column.
pool_se_fs <- function(data, se_names, reduce, group_col) {
  reduce_fn <- if (reduce == "median") stats::median else mean
  if (is.null(group_col)) {
    return(vapply(se_names, function(v) {
      reduce_fn(data[[paste0("fs_", v, "_se")]], na.rm = TRUE)
    }, numeric(1)))
  }
  gvals <- data[[group_col]]
  gnames <- fs_group_order(gvals)
  out <- do.call(rbind, lapply(gnames, function(g) {
    rows_g <- gvals == g
    vals <- vapply(se_names, function(v) {
      reduce_fn(data[[paste0("fs_", v, "_se")]][rows_g], na.rm = TRUE)
    }, numeric(1))
    data.frame(t(vals), check.names = FALSE)
  }))
  colnames(out) <- se_names
  rownames(out) <- gnames
  out
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
# (groups) in order.
tspa_schema_sf <- function(model, se) {
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
  do.call(rbind, rows)
}

# Multi-factor (fsT/fsL/fsb) schema: per latent, one struct row per score
# term and per group; error rows follow the lower triangle (incl. diagonal)
# of fsT in column-major order — the legacy per-group value routing made
# explicit and unit-testable; per-score intercept rows when fsb is given.
tspa_schema_mf <- function(model, fsT, fsL, fsb) {
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

tspa_sf <- function(model, data, se = NULL) {
  if (nrow(se) != 0) {
    return(tspa_render(tspa_schema_sf(model, se), style = "sf"))
  }
}

tspa_mf <- function(model, data, fsT, fsL, fsb) {
  tspa_render(tspa_schema_mf(model, fsT, fsL, fsb), style = "mf")
}

# ---------------------------------------------------------------------------
# Product-score (get_fs_int) auto-alias: the schema's generated model name
# for latent `v` is `fs_v`; a data column `fs_a:fs_b` (a,b latent names in
# se_fs) whose names concatenate to `v` is copied into the working data as
# `fs_v`. Old user models that pre-rename the product-score column keep
# working because the alias is a no-op when `fs_v` already exists.
# ---------------------------------------------------------------------------

tspa_sf_alias <- function(data, se) {
  is_lst <- inherits(data, "list") && !is.data.frame(data)
  dnames <- if (is_lst) names(data[[1]]) else names(data)
  se_names <- colnames(se)
  aliases <- character()
  for (v in se_names) {
    tgt <- paste0("fs_", v)
    if (tgt %in% dnames) next
    cand <- character()
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
      if (paste0(a, b) == v || paste0(b, a) == v) cand <- c(cand, col)
    }
    if (length(cand) == 0) next
    if (length(cand) > 1) {
      stop(
        "Cannot determine which product-score column in the input data ",
        "corresponds to the latent variable '", v, "': ",
        paste0("\"", cand, "\"", collapse = " or "),
        ". Rename it to \"", tgt, "\" to disambiguate."
      )
    }
    if (is_lst) {
      for (i in seq_along(data)) data[[i]][[tgt]] <- data[[i]][[cand]]
    } else {
      data[[tgt]] <- data[[cand]]
    }
    aliases <- c(aliases, paste0(tgt, " <- ", cand))
  }
  list(data = data, aliases = aliases)
}
