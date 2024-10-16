########################### Test get_fs_int function ##############################

test_that("get_fs_int() works with toy data set", code = {
  # Artificial data set
  set.seed(1230)
  dd <- data.frame(
    "a" = rnorm(5),
    "b" = runif(5),
    "w" = 1:5,
    "se_a" = rep(.15, 5),
    "se_b" = rgamma(5, 1, 1),
    "se_w" = rep(.2, 5),
    "ld_a" = 1,
    "ld_b" = c(1, 1, 1, .8, .9),
    "ld_w" = 0.7
  )

  fsint1 <- get_fs_int(dd,
    fs_name = c("a", "b", "w"), se = c("se_a", "se_b", "se_w"),
    loading = c("ld_a", "ld_b", "ld_w"), model = "a:w + a:b"
  )
  expect_equal(cor(fsint1[["a:b"]], dd[["a"]] * dd[["b"]]), 1)
  expect_in(c("a:w", "a:b", "a:b_ld", "a:w_se"), names(fsint1))
})

# Simulate a dat
sample_dat <- data.frame(
  fs1 = rnorm(100),
  fs1_se = rnorm(100),
  fs1_ld = rnorm(100),
  fs2 = rnorm(100),
  fs2_se = rnorm(100),
  fs2_ld = rnorm(100),
  fs3 = rnorm(100),
  fs3_se = rnorm(100),
  fs3_ld = rnorm(100)
)

# Test that data input must be a data frame
test_that("Data input must be a data frame", {
  expect_error(get_fs_int(dat = list(),
                          fs_name = c("fs1", "fs2", "fs3"),
                          se_fs = c("fs1_se", "fs2_se", "fs3_se"),
                          loading_fs = c("fs1_ld", "fs2_ld", "fs3_ld")))
})

# Test for correct input checks
test_that("Inputs are checked for correct types and existence", {
  expect_error(get_fs_int(dat = sample_dat,
                          fs_name = c("fs1", "abc", "fs3"),
                          se_fs = c("fs1_se", "fs2_se", "fs3_se"),
                          loading_fs = c("fs1_ld", "fs2_ld", "fs3_ld")))
})

# Test for handling numeric lat_var and its length
test_that("lat_var must be numeric and match the length of fs_name", {
  expect_error(get_fs_int(dat = sample_dat,
                          fs_name = c("fs1", "fs2", "fs3"),
                          se_fs = c("fs1_se", "fs2_se", "fs3_se"),
                          loading_fs = c("fs1_ld", "fs2_ld", "fs3_ld"),
                          lat_var = "1"))
  expect_error(get_fs_int(dat = sample_dat,
                          fs_name = c("fs1", "fs2"),
                          se_fs = c("fs1_se", "fs2_se"),
                          loading_fs = c("fs1_ld", "fs2_ld"),
                          lat_var = c(1, 2, 3)))
})

# Test function with no `lat_var` or `model`
test_that("Product indicators are correctly calculated", {
  result <- get_fs_int(dat = sample_dat,
                       fs_name = c("fs1", "fs2", "fs3"),
                       se_fs = c("fs1_se", "fs2_se", "fs3_se"),
                       loading_fs = c("fs1_ld", "fs2_ld", "fs3_ld"))
  expect_true("fs1:fs2" %in% names(result))
  expect_true("fs1:fs2_se" %in% names(result))
  expect_true("fs1:fs2_ld" %in% names(result))
  expect_true("fs1:fs3" %in% names(result))
  expect_true("fs1:fs3_se" %in% names(result))
  expect_true("fs1:fs3_ld" %in% names(result))
  expect_true("fs2:fs3" %in% names(result))
  expect_true("fs2:fs3_se" %in% names(result))
  expect_true("fs2:fs3_ld" %in% names(result))
})

