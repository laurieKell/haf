test_that("rfbParams function works correctly", {
  # Test with valid parameters
  linf <- 35
  lc <- 25
  k <- 0.25
  
  # Create example index
  set.seed(123)
  index <- FLQuant(rlnorm(20, meanlog = 10, sdlog = 0.3),
                   dimnames = list(year = 2000:2019))
  
  cntrl <- rfbParams(linf = linf, lc = lc, k = k, index = index)
  
  # Check structure
  expect_true(inherits(cntrl, "FLPar"))
  expect_true(all(c("L[F]=M", "Itrig", "multiplier") %in% dimnames(cntrl)$params))
  
  # Check calculations
  expected_lfm <- 0.75 * lc + 0.25 * linf
  expect_equal(c(cntrl["L[F]=M"]), expected_lfm, tolerance = 1e-6)
  
  expected_itrig <- 1.4 * min(index, na.rm = TRUE)
  expect_equal(c(cntrl["Itrig"]), expected_itrig, tolerance = 1e-6)
  
  # Check multiplier based on k value
  if (k < 0.2) {
    expect_equal(c(cntrl["multiplier"]), 0.95, tolerance = 1e-6)
  } else {
    expect_equal(c(cntrl["multiplier"]), 0.90, tolerance = 1e-6)
  }
})

test_that("rfbParams handles different k values", {
  linf <- 35
  lc <- 25
  
  # Create example index
  set.seed(123)
  index <- FLQuant(rlnorm(20, meanlog = 10, sdlog = 0.3),
                   dimnames = list(year = 2000:2019))
  
  # Test with k < 0.2
  cntrl_low_k <- rfbParams(linf = linf, lc = lc, k = 0.15, index = index)
  expect_equal(c(cntrl_low_k["multiplier"]), 0.95, tolerance = 1e-6)
  
  # Test with k >= 0.2
  cntrl_high_k <- rfbParams(linf = linf, lc = lc, k = 0.25, index = index)
  expect_equal(c(cntrl_high_k["multiplier"]), 0.90, tolerance = 1e-6)
})

test_that("rfbRule parameter validation works", {
  # Test input validation
  expect_error(rfbRule(iYr = "not_numeric"), "iYr must be numeric")
  expect_error(rfbRule(iYr = 2024), "missing")
}) 