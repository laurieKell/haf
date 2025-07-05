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
#' @param fCV Coefficient of variation for FMSY uncertainty (default: 0.2)
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
#' eql=eql(ple4, model="bevholtSV")
#' 
#' # Simple FMSY projection
#' fmsy_proj=Fmsy(ple4, eql, nits=50)
#' 
#' # FMSY projection with TAC bounds (20% stability clause)
#' fmsy_bounded=Fmsy(ple4, eql, nits=50, bnd=c(0.8, 1.2))
#' 
#' # Plot results
#' plot(FLStocks("FMSY"=fmsy_proj, "FMSY_Bounded"=fmsy_bounded))
#' }
#'
#' @export
setGeneric("Fmsy", function(object, eql, ...) standardGeneric("Fmsy"))

#' @rdname Fmsy
#' @export
setMethod("Fmsy", signature(object="FLStock", eql="FLBRP"),
          function(object, eql, nits=100, fmsy=NULL, 
                   maxYear=dims(object)$maxyear + 25, 
                   fCV=0.2, bnd=NULL, seed=8778, ...) {
            
            # Set seed for reproducibility
            set.seed(seed)
            
            # Get current year and extend stock to projection period
            currentYear=dims(object)$maxyear
            stk=fwdWindow(object, end=maxYear)
            
            # Extract FMSY from reference points if not provided
            if (is.null(fmsy)) {
              fmsy=refpts(eql)["msy", "harvest", drop=TRUE]
            }
            
            # Check if FMSY is available
            if (is.na(fmsy)) {
              warning("FMSY not available in reference points. Returning NULL.")
              return(NULL)
            }
            
            # Extract stock-recruitment relationship
            sr=attributes(eql)$sr
            if (is.null(sr)) {
              stop("Stock-recruitment relationship not found in equilibrium object.")
            }
            
            # Calculate recruitment deviations (following FmsyV1 exactly)
            recHat=predict(sr)
            recObs=attributes(eql)$rec.obs
            recObs=recObs[, dimnames(recObs)$year %in% dimnames(recHat)$year]
            recDevs=recObs / recHat
            recDevs=recDevs - median(recDevs)
            
            # Create recruitment deviations for projection period
            recDevs=FLQuant(
              sample(c(recDevs), dim(recDevs)[2] * nits, TRUE),
              dimnames=list(year=currentYear:maxYear, iter=seq(nits))
            )
            
            # Create FMSY with uncertainty (following FmsyV1 exactly)
            fmsy_quant=FLQuant(c(fmsy), dimnames=list(year=currentYear:maxYear))
            fmsy_quant=rlnorm(nits, log(fmsy_quant), fCV)
            
            # Perform stock projection (following FmsyV1 exactly)
            if (is.null(bnd)) {
              # Simple projection without TAC bounds
              stk=fwd(stk, f=fmsy_quant, sr=eql, deviances=exp(recDevs))
            } else {
              # Projection with TAC bounds for stability
              for (iYr in ac(currentYear:maxYear)) {
                stk=fwd(stk, f=fmsy_quant[, iYr], sr=eql, deviances=exp(recDevs))
                ctc=qmin(qmax(catch(stk)[, iYr], catch(stk)[, ac(an(iYr)-1)] * bnd[1]), 
                        catch(stk)[, ac(an(iYr)-1)] * bnd[2])
                stk=fwd(stk, catch=ctc, sr=eql, deviances=exp(recDevs))
              }
            }
            
            return(stk)
          })

#' @rdname Fmsy
#' @export
setMethod("Fmsy", signature(object="FLStocks", eql="FLBRPs"),
          function(object, eql, nits=100, fmsy=NULL, 
                   maxYear=NULL, fCV=0.2, bnd=NULL, seed=8778, ...) {
            
            # Set seed for reproducibility
            set.seed(seed)
            
            # Get common names
            stock_names=names(object)
            eql_names=names(eql)
            
            if (!all(stock_names %in% eql_names)) {
              stop("Stock names must match equilibrium object names")
            }
            
            # Set default maxYear if not provided
            if (is.null(maxYear)) {
              maxYear=max(sapply(object, function(x) dims(x)$maxyear)) + 25
            }
            
            # Perform FMSY projections for each stock
            results=lapply(stock_names, function(stock_name) {
              tryCatch({
                Fmsy(object[[stock_name]], eql[[stock_name]], 
                     nits=nits, fmsy=fmsy, maxYear=maxYear, 
                     fCV=fCV, bnd=bnd, seed=seed, ...)
              }, error=function(e) {
                warning(paste("Error projecting stock", stock_name, ":", e$message))
                return(NULL)
              })
            })
            
            # Remove NULL results and create FLStocks object
            valid_results=results[!sapply(results, is.null)]
            valid_names=stock_names[!sapply(results, is.null)]
            
            if (length(valid_results) == 0) {
              stop("No successful projections")
            }
            
            names(valid_results)=valid_names
            return(FLStocks(valid_results))
          }) 



FmsyV1<-function(stk,eql,nits=100,fmsy=benchmark(stk)["fmsy"],maxYear=dims(stk)$maxyear+25,f.cv=0.2,bnd=NA){
  set.seed(8778)
  
  currentYear=dims(stk)$maxyear
  stk    =fwdWindow(stk,end=maxYear)
  
  if (is.na(fmsy)) return(NULL)
  
  recHat =predict(attributes(eql)$sr)
  recObs =attributes(eql)$rec.obs
  recObs =recObs[,dimnames(recObs)$year%in%dimnames(recHat)$year]
  recDevs=recObs%/%recHat  
  recDevs=recDevs-median(recDevs)
  
  recDevs=FLQuant(sample(c(recDevs),dim(recDevs)[2]*nits,TRUE),
                  dimnames=list(year=currentYear:maxYear,iter=seq(nits)))
  
  fmsy   =FLQuant(c(fmsy),dimnames=list(year=currentYear:maxYear))
  fmsy   =rlnorm(nits,log(fmsy),f.cv)
  
  if (any(is.na(bnd)))
    stk=fwd(stk,f=fmsy,sr=eql,deviances=exp(recDevs))
  else{
    for (iYr in ac((currentYear):maxYear)){
      stk=fwd(stk,f=fmsy[,iYr],sr=eql,deviances=exp(recDevs))
      ctc=qmin(qmax(catch(stk)[,iYr],catch(stk)[,ac(an(iYr)-1)]*bnd[1]),catch(stk)[,ac(an(iYr)-1)]*bnd[2])
      stk=fwd(stk,catch=ctc,sr=eql,deviances=exp(recDevs))}}
  
  stk}  

# Test function to debug Fmsy issues
testFmsy <- function(stk, eql, nits=50) {
  cat("Testing Fmsy method...\n")
  cat("Stock dimensions:", dim(stk), "\n")
  cat("Stock years:", range(dimnames(stk)$year), "\n")
  cat("FMSY value:", refpts(eql)["msy", "harvest"], "\n")
  
  # Try the working FmsyV1 function first
  cat("Testing FmsyV1...\n")
  result_v1 <- try(FmsyV1(stk, eql, nits=nits))
  if (!inherits(result_v1, "try-error")) {
    cat("FmsyV1 worked successfully!\n")
  } else {
    cat("FmsyV1 failed:", result_v1, "\n")
  }
  
  # Try the new Fmsy method
  cat("Testing Fmsy method...\n")
  result_new <- try(Fmsy(stk, eql, nits=nits))
  if (!inherits(result_new, "try-error")) {
    cat("Fmsy method worked successfully!\n")
  } else {
    cat("Fmsy method failed:", result_new, "\n")
  }
  
  return(list(v1=result_v1, new=result_new))
}  