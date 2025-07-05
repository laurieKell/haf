#' Management Strategy Evaluation shortcut Method
#'
#' @description
#' Performs Management Strategy Evaluation (MSE) using Monte Carlo simulations
#' to assess the performance of harvest control rules under uncertainty.
#' This method implements a "shortcut" MSE approach that uses historical
#' recruitment residuals and observation errors to simulate future scenarios.
#'
#' @param object An FLStock or FLStocks object
#' @param eql An FLBRP or FLBRPs object with fitted stock-recruitment relationships
#' @param endYear Final year for projections (default: 2050)
#' @param nits Number of Monte Carlo iterations (default: 100)
#' @param seed Random seed for reproducibility (default: 1233)
#' @param hcrParams HCR parameters (if NULL, uses benchmark values)
#' @param bndTac TAC bounds for stability constraints (default: c(0.8, 1.25))
#' @param workers Number of workers for parallel processing (default: availableCores()-2)
#' @param ... Additional arguments
#'
#' @return A list containing MSE results for each stock
#'
#' @details
#' The shortcut method performs the following steps:
#' 
#' 1. **Data Preparation**: Ensures positive values for catch, discards, landings, and stock numbers
#' 2. **Recruitment Analysis**: Extracts and samples historical recruitment residuals
#' 3. **Observation Error**: Generates observation errors for future projections
#' 4. **HCR Parameters**: Sets up harvest control rule parameters from benchmark values
#' 5. **Parallel Simulation**: Runs Monte Carlo simulations across multiple workers
#' 6. **Results Collection**: Compiles and returns MSE results
#'
#' The method supports both single stocks and multiple stocks, and can use
#' different stock-recruitment relationships (Beverton-Holt, Ricker, etc.).
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ple4)
#' 
#' # Fit stock-recruitment relationship
#' eql <- eql(ple4, model = "bevholtSV")
#' 
#' # Run MSE for single stock
#' mse_result <- MSE(ple4, eql, nits = 50, endYear = 2030)
#' 
#' # Run MSE for multiple stocks
#' stocks <- FLStocks("stock1" = ple4, "stock2" = ple4)
#' eqls <- FLBRPs("stock1" = eql, "stock2" = eql)
#' mse_results <- MSE(stocks, eqls, nits = 50)
#' }
#'
#' @export
setGeneric("shortcut", function(object, eql, ...) standardGeneric("shortcut"))

#' @rdname shortcut
#' @export
setMethod("shortcut", signature(object = "FLStock", eql = "FLBRP"),
          function(object, eql, endYear = 2050, nits = 100, seed = 1233,
                   hcrParams = NULL, bndTac = c(0.8, 1.25), workers = NULL, ...) {
            
            # Set up parallel processing with future.apply
            if (is.null(workers)) {
              workers <- max(1, future::availableCores() - 2)
            }
            
            # Set up future plan
            old_plan <- future::plan(future::multisession, workers = workers)
            on.exit(future::plan(old_plan), add = TRUE)
            
            # Set seed for reproducibility
            set.seed(seed)
            
            # Ensure positive values for stock data
            object <- ensurePositive(object)
            
            # Check if FMSY is available
            fmsy <- benchmark(object, eql)["fmsy"]
            if (is.na(fmsy)) {
              warning("FMSY not available in benchmark. Returning NULL.")
              return(NULL)
            }
            
            # Extract recruitment residuals
            recResiduals <- attributes(eql)$rec.residuals
            if (is.null(recResiduals)) {
              warning("Recruitment residuals not found in equilibrium object.")
              return(NULL)
            }
            
            # Use last 20 years of recruitment residuals
            yrs <- recResiduals[, rev(dim(recResiduals)[2] - 0:min(dim(recResiduals)[2] - 1, 20))]
            recDevs <- recResiduals[, rev(dim(yrs)[2] - 0:min(dim(recResiduals)[2] - 1, 20))]
            recDevs <- recDevs - median(recDevs)
            
            # Create recruitment deviations for projection period
            recDevs <- FLQuant(
              sample(c(recDevs), dim(recDevs)[2] * nits, TRUE),
              dimnames = list(year = dims(object)$maxyear:endYear, iter = seq(nits))
            )
            
            # Generate observation errors (shortcut approach)
            # Note: This uses a simplified approach - in practice, you might want more sophisticated error generation
            obsError <- rlnorm(nits, 
                              FLQuant(0, dimnames = list(year = 2019:endYear)), 
                              sqrt(var(recDevs)))
            
            # Set up HCR parameters if not provided
            if (is.null(hcrParams)) {
              hcrParams <- FLPar(
                ftar = c(benchmark(object, eql)["fmsy"]),
                btrig = c(benchmark(object, eql)["btrigger"]),
                fmin = 0.005,
                bmin = 0,
                blim = benchmark(object, eql)["blim"]
              )
            }
            
            # Prepare stock for projection
            stk <- propagate(fwdWindow(object, end = endYear + 1), nits)
            
            # Run MSE simulation
            result <- tryCatch({
              hcrICES(stk, eql, exp(recDevs), hcrParams,
                     start = dims(object)$maxyear, end = endYear,
                     interval = 1, lag = 1, err = obsError, bndTac = bndTac)
            }, error = function(e) {
              warning(paste("MSE simulation failed:", e$message))
              return(NULL)
            })
            
            return(result)
          })

#' @rdname shortcut
#' @export
setMethod("shortcut", signature(object = "FLStocks", eql = "FLBRPs"),
          function(object, eql, endYear = 2050, nits = 100, seed = 1233,
                   hcrParams = NULL, bndTac = c(0.8, 1.25), workers = NULL, ...) {
            
            # Set up parallel processing with future.apply
            if (is.null(workers)) {
              workers <- max(1, future::availableCores() - 2)
            }
            
            # Set up future plan
            old_plan <- future::plan(future::multisession, workers = workers)
            on.exit(future::plan(old_plan), add = TRUE)
            
            # Get stock names
            stockNames <- names(object)
            eqlNames <- names(eql)
            
            if (!all(stockNames %in% eqlNames)) {
              stop("Stock names must match equilibrium object names")
            }
            
            # Set seed for reproducibility
            set.seed(seed)
            
            # Create MSE function for single stock
            runSingleMSE <- function(id) {
              # Ensure positive values for stock data
              stk <- ensurePositive(object[[id]])
              
              # Check if FMSY is available
              fmsy <- benchmark(stk, eql[[id]])["fmsy"]
              if (is.na(fmsy)) return(NULL)
              
              # Extract recruitment residuals
              recResiduals <- attributes(eql[[id]])$rec.residuals
              if (is.null(recResiduals)) return(NULL)
              
              # Use last 20 years of recruitment residuals
              yrs <- recResiduals[, rev(dim(recResiduals)[2] - 0:min(dim(recResiduals)[2] - 1, 20))]
              recDevs <- recResiduals[, rev(dim(yrs)[2] - 0:min(dim(recResiduals)[2] - 1, 20))]
              recDevs <- recDevs - median(recDevs)
              
              # Create recruitment deviations for projection period
              recDevs <- FLQuant(
                sample(c(recDevs), dim(recDevs)[2] * nits, TRUE),
                dimnames = list(year = dims(stk)$maxyear:endYear, iter = seq(nits))
              )
              
              # Generate observation errors
              obsError <- rlnorm(nits, 
                                FLQuant(0, dimnames = list(year = 2019:endYear)), 
                                sqrt(var(recDevs)))
              
              # Set up HCR parameters if not provided
              if (is.null(hcrParams)) {
                hcrParams <- FLPar(
                  ftar = c(benchmark(stk, eql[[id]])["fmsy"]),
                  btrig = c(benchmark(stk, eql[[id]])["btrigger"]),
                  fmin = 0.005,
                  bmin = 0,
                  blim = benchmark(stk, eql[[id]])["blim"]
                )
              }
              
              # Prepare stock for projection
              stk <- propagate(fwdWindow(stk, end = endYear + 1), nits)
              
              # Run MSE simulation
              tryCatch({
                hcrICES(stk, eql[[id]], exp(recDevs), hcrParams,
                       start = dims(stk)$maxyear, end = endYear,
                       interval = 1, lag = 1, err = obsError, bndTac = bndTac)
              }, error = function(e) {
                warning(paste("MSE simulation failed for stock", id, ":", e$message))
                return(NULL)
              })
            }
            
            # Run MSE for each stock in parallel using future.apply
            results <- future.apply::future_lapply(stockNames, runSingleMSE)
            
            # Name results
            names(results) <- stockNames
            
            # Remove NULL results
            validResults <- results[!sapply(results, is.null)]
            
            if (length(validResults) == 0) {
              stop("No successful MSE simulations")
            }
            
            return(validResults)
          })

#' Helper function to ensure positive values in FLStock
#' @param stk FLStock object
#' @return FLStock with positive values
ensurePositive <- function(stk) {
  # Ensure positive values for various components
  catch.n(stk)[catch.n(stk) <= 0] <- 0.0001
  discards.n(stk)[discards.n(stk) <= 0] <- 0.0001
  landings.n(stk)[landings.n(stk) <= 0] <- 0.0001
  stock.n(stk)[stock.n(stk) <= 0] <- 0.0001
  
  return(stk)
}

#' Benchmark function to extract reference points
#' @param stk FLStock object
#' @param eql FLBRP object
#' @return Named vector with benchmark values
benchmark <- function(stk, eql) {
  # Extract reference points from equilibrium object
  refs <- refpts(eql)
  
  # Create benchmark vector with common reference points
  bench <- c(
    fmsy = refs["msy", "harvest", drop = TRUE],
    btrigger = refs["msy", "ssb", drop = TRUE] * 0.8,  # Common default
    blim = refs["msy", "ssb", drop = TRUE] * 0.5       # Common default
  )
  
  return(bench)
}

