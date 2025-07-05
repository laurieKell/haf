#' Management Strategy Evaluation (MSE) Method
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
#' @param cores Number of CPU cores for parallel processing (default: detectCores()-4)
#' @param ... Additional arguments
#'
#' @return A list containing MSE results for each stock
#'
#' @details
#' The MSE method performs the following steps:
#' 
#' 1. **Data Preparation**: Ensures positive values for catch, discards, landings, and stock numbers
#' 2. **Recruitment Analysis**: Extracts and samples historical recruitment residuals
#' 3. **Observation Error**: Generates observation errors for future projections
#' 4. **HCR Parameters**: Sets up harvest control rule parameters from benchmark values
#' 5. **Parallel Simulation**: Runs Monte Carlo simulations across multiple cores
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
setGeneric("MSE", function(object, eql, ...) standardGeneric("MSE"))

#' @rdname MSE
#' @export
setMethod("MSE", signature(object = "FLStock", eql = "FLBRP"),
          function(object, eql, endYear = 2050, nits = 100, seed = 1233,
                   hcrParams = NULL, bndTac = c(0.8, 1.25), cores = NULL, ...) {
            
            # Set up parallel processing
            if (is.null(cores)) {
              cores <- max(1, parallel::detectCores() - 4)
            }
            
            # Register parallel backend
            cl <- parallel::makeCluster(cores)
            doParallel::registerDoParallel(cl)
            
            # Set seed for reproducibility
            set.seed(seed)
            
            # Ensure positive values for stock data
            object <- ensurePositive(object)
            
            # Check if FMSY is available
            fmsy <- benchmark(object)["fmsy"]
            if (is.na(fmsy)) {
              warning("FMSY not available in benchmark. Returning NULL.")
              parallel::stopCluster(cl)
              return(NULL)
            }
            
            # Extract recruitment residuals
            recResiduals <- attributes(eql)$rec.residuals
            if (is.null(recResiduals)) {
              warning("Recruitment residuals not found in equilibrium object.")
              parallel::stopCluster(cl)
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
                ftar = c(benchmark(object)["fmsy"]),
                btrig = c(benchmark(object)["btrigger"]),
                fmin = 0.005,
                bmin = 0,
                blim = benchmark(object)["blim"]
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
            
            # Clean up parallel processing
            parallel::stopCluster(cl)
            
            return(result)
          })

#' @rdname MSE
#' @export
setMethod("MSE", signature(object = "FLStocks", eql = "FLBRPs"),
          function(object, eql, endYear = 2050, nits = 100, seed = 1233,
                   hcrParams = NULL, bndTac = c(0.8, 1.25), cores = NULL, ...) {
            
            # Set up parallel processing
            if (is.null(cores)) {
              cores <- max(1, parallel::detectCores() - 4)
            }
            
            # Register parallel backend
            cl <- parallel::makeCluster(cores)
            doParallel::registerDoParallel(cl)
            
            # Get stock names
            stockNames <- names(object)
            eqlNames <- names(eql)
            
            if (!all(stockNames %in% eqlNames)) {
              stop("Stock names must match equilibrium object names")
            }
            
            # Set seed for reproducibility
            set.seed(seed)
            
            # Run MSE for each stock in parallel
            results <- foreach(id = stockNames,
                             .combine = list,
                             .maxcombine = length(object),
                             .packages = c("FLCore", "FLBRP", "FLasher", "FLCandy", "haf")
            ) %dopar% {
              
              # Ensure positive values for stock data
              stk <- ensurePositive(object[[id]])
              
              # Check if FMSY is available
              fmsy <- benchmark(stk)["fmsy"]
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
                  ftar = c(benchmark(stk)["fmsy"]),
                  btrig = c(benchmark(stk)["btrigger"]),
                  fmin = 0.005,
                  bmin = 0,
                  blim = benchmark(stk)["blim"]
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
            
            # Clean up parallel processing
            parallel::stopCluster(cl)
            
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

ncore=detectCores()-4  
registerDoParallel(ncore)

catch.n(   icesdata[[46]])[catch.n(   icesdata[[46]])<=0]=0.0001
discards.n(icesdata[[46]])[discards.n(icesdata[[46]])<=0]=0.0001
landings.n(icesdata[[46]])[landings.n(icesdata[[46]])<=0]=0.0001
stock.n(   icesdata[[46]])[stock.n(   icesdata[[46]])<=0]=0.0001

mseBH=foreach(id=stks[!stks[,1]%in%BH,".id"],#names(icesdata),
                     .combine   =list,
                     .maxcombine=length(icesdata),
                     .packages  =c("FLCore","FLBRP","FLasher","FLCandy")
                     ) %dopar% {
    source("C:/flrpapers/erp/source/hcrICESV3.R")

    endYear=2050
    nits   =100
    set.seed(1233)

    shortCut=rlnorm(nits,FLQuant(0,dimnames=list(year=2019:endYear)),var(err[-1,1])^0.5)
      
    stk=icesdata[[id]]
      
    if (is.na(c(benchmark(stk)["fmsy"]))) return(NULL)
     
    yrs=attributes(bhs[[id]])$rec.residuals
    yrs=yrs[,rev(dim(yrs)[2]-0:min(dim(yrs)[2]-1,20))]

    recDevs=attributes(bhs[[id]])$rec.residuals
    recDevs=recDevs[,rev(dim(yrs)[2]-0:min(dim(recDevs)[2]-1,20))]
    recDevs=recDevs-median(recDevs)
    recDevs=FLQuant(sample(c(recDevs),dim(recDevs)[2]*nits,TRUE),
                      dimnames=list(year=dims(icesdata[[id]])$maxyear:endYear,iter=seq(nits)))

      
    hcrPar=FLPar(c(ftar =c(benchmark(stk)["fmsy"]),
                   btrig=c(benchmark(stk)["btrigger"]),
                   fmin =FLPar(fmin=0.005),
                   bmin =FLPar(bmin=0),
                   blim =benchmark(stk)["blim"]))
      
    stk=propagate(fwdWindow(stk,end=endYear+1),nits)
    rtn=tryIt(hcrICES(stk,bhs[[id]],exp(recDevs),
                hcrPar,
                start=dims(icesdata[[id]])$maxyear,end=endYear,
                interval=1,lag=1,
                err=shortCut,
                bndTac=c(0.8,1.25)))
    
    save(rtn,file=paste(file.path("../data/results",id),".BH.RData"))
    
    rtn}

mseRK=foreach(id=stks[!stks[,1]%in%RK,".id"],#names(icesdata),
                     .combine   =list,
                     .maxcombine=length(icesdata),
                     .packages  =c("FLCore","FLBRP","FLasher","FLCandy")
                     ) %dopar% {
      
    source("C:/flrpapers/erp/source/hcrICESV3.R")
    endYear=2050
    nits   =100
    set.seed(1233)

    shortCut=rlnorm(nits,FLQuant(0,dimnames=list(year=2019:endYear)),var(err[-1,1])^0.5)
      
    stk=icesdata[[id]]
      
    if (is.na(c(benchmark(stk)["fmsy"]))) return(NULL)
      
    recDevs=attributes(rks[[id]])$rec.residuals
    recDevs=recDevs[,rev(dim(yrs)[2]-0:min(dim(recDevs)[2]-1,20))]
    recDevs=recDevs-median(recDevs)
    recDevs=FLQuant(sample(c(recDevs),dim(recDevs)[2]*nits,TRUE),
                      dimnames=list(year=dims(icesdata[[id]])$maxyear:endYear,iter=seq(nits)))
          
    hcrPar=FLPar(c(ftar =c(benchmark(stk)["fmsy"]),
                   btrig=c(benchmark(stk)["btrigger"]),
                   fmin =FLPar(fmin=0.005),
                   bmin =FLPar(bmin=0),
                   blim =benchmark(stk)["blim"]))
      
    stk=propagate(fwdWindow(stk,end=endYear+1),nits)
    rtn=tryIt(hcrICES(stk,rks[[id]],exp(recDevs),
                hcrPar,
                start=dims(icesdata[[id]])$maxyear,end=endYear,
                interval=1,lag=1,
                err=shortCut,
                bndTac=c(0.8,1.25)))
    
    save(rtn,file=paste(file.path("../data/results",id),".RK.RData"))
    
    rtn}

stopImplicitCluster()

