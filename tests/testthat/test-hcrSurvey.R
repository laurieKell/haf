test_that("hcrSurvey function works correctly", {
  # Test with valid parameters
  iYr <- 2024
  lag <- 1
  
  # Create example data
  index <- FLQuant(seq(1000, 8000, length.out = 25),
                   dimnames = list(year = 2000:2024))
  catch <- FLQuant(800, dimnames = list(year = 2000:2024))
  cntrl <- FLPar(fmsy = 0.16, btrig = 5000, bbuf = 2500)
  
  # Test function execution
  result <- hcrSurvey(iYr = iYr, index = index, catch = catch, 
                      cntrl = cntrl, lag = lag)
  
  # Check structure
  expect_true(inherits(result, "FLQuant"))
  expect_equal(dimnames(result)$year, as.character(iYr))
})

test_that("hcrSurvey handles different zones correctly", {
  iYr <- 2024
  lag <- 1
  catch <- FLQuant(800, dimnames = list(year = 2000:2024))
  cntrl <- FLPar(fmsy = 0.16, btrig = 5000, bbuf = 2500)
  
  # Test green zone (B >= Btrigger)
  index_green <- FLQuant(6000, dimnames = list(year = 2000:2024))
  result_green <- hcrSurvey(iYr = iYr, index = index_green, catch = catch, 
                            cntrl = cntrl, lag = lag)
  
  # Test yellow zone (Bbuffer <= B < Btrigger)
  index_yellow <- FLQuant(3500, dimnames = list(year = 2000:2024))
  result_yellow <- hcrSurvey(iYr = iYr, index = index_yellow, catch = catch, 
                             cntrl = cntrl, lag = lag)
  
  # Test red zone (B < Bbuffer)
  index_red <- FLQuant(1500, dimnames = list(year = 2000:2024))
  result_red <- hcrSurvey(iYr = iYr, index = index_red, catch = catch, 
                          cntrl = cntrl, lag = lag)
  
  # All results should be FLQuant objects
  expect_true(inherits(result_green, "FLQuant"))
  expect_true(inherits(result_yellow, "FLQuant"))
  expect_true(inherits(result_red, "FLQuant"))
})

test_that("hcrSurvey error handling works", {
  iYr <- 2024
  lag <- 1
  index <- FLQuant(1000, dimnames = list(year = 2000:2024))
  catch <- FLQuant(800, dimnames = list(year = 2000:2024))
  
  # Test with invalid parameters (negative values)
  cntrl_invalid <- FLPar(fmsy = 0.16, btrig = -5000, bbuf = 2500)
  expect_error(hcrSurvey(iYr = iYr, index = index, catch = catch, 
                         cntrl = cntrl_invalid, lag = lag),
               "Biomass values must be positive")
  
  # Test with negative index
  index_negative <- FLQuant(-1000, dimnames = list(year = 2000:2024))
  cntrl_valid <- FLPar(fmsy = 0.16, btrig = 5000, bbuf = 2500)
  expect_error(hcrSurvey(iYr = iYr, index = index_negative, catch = catch, 
                         cntrl = cntrl_valid, lag = lag),
               "Biomass values must be positive")
}) 