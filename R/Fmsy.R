#' FMSY Projection Method
#'
#' @description
#' Performs FMSY (Maximum Sustainable Yield fishing mortality) projections for FLR stock objects.
#' This method projects stock dynamics under FMSY fishing mortality with uncertainty analysis
#' through Monte Carlo simulations.
#'
#' @param object An FLStock object
#' @param eql An FLBRP object with fitted stock-recruitment relationship
#' @param nits Number of iterations for uncertainty analysis (default: 100)
#' @param fmsy FMSY value (if NULL, uses reference points from eql)
#' @param maxYear Maximum projection year (default: current year + 25)
#' @param f_cv Coefficient of variation for FMSY uncertainty (default: 0.2)
#' @param bnd TAC bounds c(lower, upper) for stability constraints (default: NULL)
#' @param seed Random seed for reproducibility (default: 8778)
#' @param ... Additional arguments
#'
#' @return An FLStock object with FMSY projections
#'
#' @details
#' The FMSY projection method performs the following steps:
#' 
#' 1. **Extract FMSY**: Uses reference points from the equilibrium object or provided value
#' 2. **Recruitment Analysis**: Calculates recruitment deviations from fitted SRR
#' 3. **Uncertainty Sampling**: Generates Monte Carlo samples for FMSY and recruitment
#' 4. **Stock Projection**: Projects stock dynamics using FLR's forward simulation
#' 5. **Stability Constraints**: Optionally applies TAC bounds to prevent unrealistic changes
#'
#' The method handles both simple projections (without TAC bounds) and constrained
#' projections with stability clauses to simulate real-world management constraints.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ple4)
#' 
#' # Fit stock-recruitment relationship
#' eql <- eql(ple4, model = "bevholtSV")
#' 
#' # Simple FMSY projection
#' fmsy_proj <- Fmsy(ple4, eql, nits = 50)
#' 
#' # FMSY projection with TAC bounds (20% stability clause)
#' fmsy_bounded <- Fmsy(ple4, eql, nits = 50, bnd = c(0.8, 1.2))
#' 
#' # Plot results
#' plot(FLStocks("FMSY" = fmsy_proj, "FMSY_Bounded" = fmsy_bounded))
#' }
#'
#' @export
setGeneric("Fmsy", function(object, eql, ...) standardGeneric("Fmsy"))

#' @rdname Fmsy
#' @export
setMethod("Fmsy", signature(object = "FLStock", eql = "FLBRP"),
          function(object, eql, nits = 100, fmsy = NULL, 
                   maxYear = dims(object)$maxyear + 25, 
                   f_cv = 0.2, bnd = NULL, seed = 8778, ...) {
            
            # Set seed for reproducibility
            set.seed(seed)
            
            # Get current year and extend stock to projection period
            currentYear <- dims(object)$maxyear
            stk <- window(object, end = maxYear)
            
            # Extract FMSY from reference points if not provided
            if (is.null(fmsy)) {
              fmsy <- refpts(eql)["msy", "harvest", drop = TRUE]
            }
            
            # Check if FMSY is available
            if (is.na(fmsy)) {
              warning("FMSY not available in reference points. Returning NULL.")
              return(NULL)
            }
            
            # Extract stock-recruitment relationship
            sr_obj <- attributes(eql)$sr
            if (is.null(sr_obj)) {
              stop("Stock-recruitment relationship not found in equilibrium object.")
            }
            
            # Calculate recruitment deviations
            rec_hat <- predict(sr_obj)
            rec_obs <- attributes(eql)$rec.obs
            
            if (is.null(rec_obs)) {
              stop("Recruitment observations not found in equilibrium object.")
            }
            
            # Align recruitment data with predicted values
            rec_obs <- rec_obs[, dimnames(rec_obs)$year %in% dimnames(rec_hat)$year]
            rec_devs <- rec_obs / rec_hat
            rec_devs <- rec_devs - median(rec_devs)
            
            # Create recruitment deviations for projection period
            rec_devs <- FLQuant(
              sample(c(rec_devs), dim(rec_devs)[2] * nits, TRUE),
              dimnames = list(year = currentYear:maxYear, iter = seq(nits))
            )
            
            # Create FMSY with uncertainty
            fmsy_quant <- FLQuant(c(fmsy), dimnames = list(year = currentYear:maxYear))
            fmsy_quant <- rlnorm(nits, log(fmsy_quant), f_cv)
            
            # Ensure proper dimensions for FMSY quant
            fmsy_quant <- propagate(fmsy_quant, nits)
            
            # Ensure recruitment deviations have proper dimensions
            rec_devs <- propagate(rec_devs, nits)
            
            # Perform stock projection
            if (is.null(bnd)) {
              # Simple projection without TAC bounds
              stk <- fwd(stk, f = fmsy_quant, sr = eql, residuals = exp(rec_devs))
            } else {
              # Projection with TAC bounds for stability
              if (length(bnd) != 2 || bnd[1] >= bnd[2]) {
                stop("TAC bounds must be a vector of length 2 with lower < upper")
              }
              
              for (iYr in ac(currentYear:maxYear)) {
                # Project with FMSY
                stk <- fwd(stk, f = fmsy_quant[, iYr], sr = eql, residuals = exp(rec_devs))
                
                # Apply TAC bounds
                prev_catch <- catch(stk)[, ac(as.numeric(iYr) - 1)]
                current_catch <- catch(stk)[, iYr]
                
                bounded_catch <- qmin(
                  qmax(current_catch, prev_catch * bnd[1]), 
                  prev_catch * bnd[2]
                )
                
                # Re-project with bounded catch
                stk <- fwd(stk, catch = bounded_catch, sr = eql, residuals = exp(rec_devs))
              }
            }
            
            return(stk)
          })

#' @rdname Fmsy
#' @export
setMethod("Fmsy", signature(object = "FLStocks", eql = "FLBRPs"),
          function(object, eql, nits = 100, fmsy = NULL, 
                   maxYear = NULL, f_cv = 0.2, bnd = NULL, seed = 8778, ...) {
            
            # Set seed for reproducibility
            set.seed(seed)
            
            # Get common names
            stock_names <- names(object)
            eql_names <- names(eql)
            
            if (!all(stock_names %in% eql_names)) {
              stop("Stock names must match equilibrium object names")
            }
            
            # Set default maxYear if not provided
            if (is.null(maxYear)) {
              maxYear <- max(sapply(object, function(x) dims(x)$maxyear)) + 25
            }
            
            # Perform FMSY projections for each stock
            results <- lapply(stock_names, function(stock_name) {
              tryCatch({
                Fmsy(object[[stock_name]], eql[[stock_name]], 
                     nits = nits, fmsy = fmsy, maxYear = maxYear, 
                     f_cv = f_cv, bnd = bnd, seed = seed, ...)
              }, error = function(e) {
                warning(paste("Error projecting stock", stock_name, ":", e$message))
                return(NULL)
              })
            })
            
            # Remove NULL results and create FLStocks object
            valid_results <- results[!sapply(results, is.null)]
            valid_names <- stock_names[!sapply(results, is.null)]
            
            if (length(valid_results) == 0) {
              stop("No successful projections")
            }
            
            names(valid_results) <- valid_names
            return(FLStocks(valid_results))
          }) 