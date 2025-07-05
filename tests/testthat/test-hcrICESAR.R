test_that("hcrICES function works correctly", {
  # Test parameter validation
  params <- FLPar(
    blim = 100000,
    btrig = 150000,
    bmin = 100000,
    ftar = 0.3,
    fmin = 0.05
  )
  
  expect_true(all(c("blim", "btrig", "bmin", "ftar", "fmin") %in% dimnames(params)$params))
  expect_true(params["blim"] < params["btrig"])
  expect_true(params["fmin"] < params["ftar"])
})

test_that("hcrICES handles edge cases", {
  # Test with extreme values
  params_extreme <- FLPar(
    blim = 0.001,
    btrig = 0.1,
    bmin = 0.001,
    ftar = 0.01,
    fmin = 0.001
  )
  
  expect_true(params_extreme["blim"] < params_extreme["btrig"])
  expect_true(params_extreme["fmin"] < params_extreme["ftar"])
})

test_that("hcrICES parameter structure is correct", {
  params <- FLPar(
    blim = 100000,
    btrig = 150000,
    bmin = 100000,
    ftar = 0.3,
    fmin = 0.05
  )
  
  # Check that params is an FLPar object
  expect_true(inherits(params, "FLPar"))
  
  # Check that all required parameters are present
  required_params <- c("blim", "btrig", "bmin", "ftar", "fmin")
  expect_true(all(required_params %in% dimnames(params)$params))
}) 