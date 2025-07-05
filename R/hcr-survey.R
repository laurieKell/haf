#' Survey-Based Harvest Control Rule
#' 
#' @title Survey-Based Harvest Control Rule (HCR)
#' 
#' @description Implements a generic survey-based harvest control rule for stocks where biomass estimates come from surveys rather than analytical assessments. This rule is particularly suitable for stocks like Nephrops, crustaceans, and other species managed through survey indices.
#' 
#' @param iYr \code{numeric} Current management year
#' @param index \code{FLQuant} Survey biomass index time series
#' @param catch \code{FLQuant} Historical catch time series
#' @param cntrl \code{FLPar} Control parameters containing fmsy, btrig, and bbuf
#' @param minCatch \code{numeric} Minimum catch constraint as proportion of previous catch (default: 0.1)
#' @param lag \code{numeric} Lag between survey and advice (default: 1)
#' 
#' @return \code{FLQuant} containing new catch advice for year \code{iYr}
#' 
#' @details
#' The survey-based HCR implements a three-zone approach:
#' \itemize{
#'   \item \strong{Green zone} (B >= Btrigger): F = FMSY
#'   \item \strong{Yellow zone} (Bbuffer <= B < Btrigger): F = FMSY × reduction factor
#'   \item \strong{Red zone} (B < Bbuffer): F = 0.8 × previous catch
#' }
#' 
#' The reduction factor in the yellow zone is calculated as:
#' \deqn{reduction = \frac{B - B_{buffer}}{B_{trigger} - B_{buffer}}}
#' 
#' This approach is particularly useful for:
#' \itemize{
#'   \item Nephrops and other crustacean stocks
#'   \item Stocks with limited analytical assessments
#'   \item Species managed primarily through survey indices
#'   \item Data-limited stocks requiring precautionary management
#' }
#' 
#' @examples
#' \dontrun{
#' # Set up control parameters for a survey-based stock
#' cntrl <- FLPar(
#'   fmsy = 0.16,    # Target fishing mortality
#'   btrig = 5000,   # Biomass trigger reference point
#'   bbuf = 2500     # Biomass buffer reference point
#' )
#' 
#' # Create example index and catch data
#' index <- FLQuant(seq(1, 10000, length.out=52), 
#'                  dimnames=list(year=1990:2041))
#' catch <- index %*% cntrl["fmsy"]
#' 
#' # Project forward using HCR
#' for (i in 1991:2040) {
#'   catch[,ac(i)] <- hcrSurvey(i, index, catch, cntrl, lag=0)
#' }
#' }
#' 
#' @export
#' @importFrom FLCore ebiomass
hcrSurvey<-function(iYr,index,catch,cntrl,minCatch = 0.1,lag=1){
  
  # Error checking
  if (any(c(index[,ac(iYr-lag)], cntrl["btrig"], cntrl["bbuf"]) <= 0)) 
    stop("Biomass values must be positive")
  
  ## Decision rule
  # default
  abc=FLQuant(c(index[,ac(iYr-lag)]%*%cntrl["fmsy"]),dimnames=dimnames(catch[,ac(iYr)]))
  
  btrig=index[,ac(iYr-lag)]<c(cntrl["btrig"])
  bbuf =index[,ac(iYr-lag)]<c(cntrl["bbuf"])
  
  ## If stock below bbuffer reduce catch
  if (any(btrig))
    abc[,,,,,bbuf]=catch[,ac(iYr-lag),,,,bbuf]*0.8
  
  ## If stock below bbuffer reduce catch
  if (any(bbuf)){
    redFactor=(index[,ac(iYr-lag)]-cntrl["bbuf"])%/%(cntrl["btrig"]%-%cntrl["bbuf"])
    abc[,,,,,btrig]=index[,ac(iYr-lag),,,,btrig]%*%cntrl["fmsy"]%*%redFactor[,,,,,btrig]}
  
  # Apply minimum catch constraint
  abc=qmax(abc,catch[,ac(iYr-lag)]*minCatch)
  
  return(abc)}

#' Nephrops-Specific Harvest Control Rule
#' 
#' @title Nephrops Harvest Control Rule (HCR)
#' 
#' @description Implements the ICES harvest control rule specifically for Nephrops stocks based on survey indices. This is a wrapper around the generic survey-based HCR with Nephrops-specific parameter defaults and documentation.
#' 
#' @param iYr \code{numeric} Current management year
#' @param index \code{FLQuant} Survey biomass index time series
#' @param catch \code{FLQuant} Historical catch time series
#' @param cntrl \code{FLPar} Control parameters containing fmsy, btrig, and bbuf
#' @param minCatch \code{numeric} Minimum catch constraint as proportion of previous catch (default: 0.1)
#' @param lag \code{numeric} Lag between survey and advice (default: 1)
#' 
#' @return \code{FLQuant} containing new catch advice for year \code{iYr}
#' 
#' @details
#' This function is a wrapper around \code{\link{hcrSurvey}} specifically for Nephrops stocks.
#' The Nephrops HCR implements a three-zone approach as described in ICES guidelines:
#' \itemize{
#'   \item \strong{Green zone} (B >= Btrigger): F = FMSY
#'   \item \strong{Yellow zone} (Bbuffer <= B < Btrigger): F = FMSY × reduction factor
#'   \item \strong{Red zone} (B < Bbuffer): F = 0.8 × previous catch
#' }
#' 
#' @examples
#' \dontrun{
#' # Set up control parameters for Nephrops
#' cntrl <- FLPar(
#'   fmsy = 0.16,    # Target fishing mortality
#'   btrig = 5000,   # Biomass trigger reference point
#'   bbuf = 2500     # Biomass buffer reference point
#' )
#' 
#' # Create example index and catch data
#' index <- FLQuant(seq(1, 10000, length.out=52), 
#'                  dimnames=list(year=1990:2041))
#' catch <- index %*% cntrl["fmsy"]
#' 
#' # Project forward using HCR
#' for (i in 1991:2040) {
#'   catch[,ac(i)] <- hcrNephrops(i, index, catch, cntrl, lag=0)
#' }
#' }
#' 
#' @export
hcrNephrops <- function(iYr, index, catch, cntrl, minCatch = 0.1, lag = 1) {
  # Call the generic survey-based HCR
  hcrSurvey(iYr = iYr, index = index, catch = catch, 
            cntrl = cntrl, minCatch = minCatch, lag = lag)
}

#' Survey-Based HCR Visualization
#' 
#' @title Survey-Based HCR Visualization
#' 
#' @description Creates a diagnostic plot showing the survey-based HCR zones and catch response
#' 
#' @param cntrl \code{FLPar} Control parameters
#' @param index \code{FLQuant} Index values for plotting
#' @param catch \code{FLQuant} Catch values for plotting
#' @param title \code{character} Plot title (default: "Survey-Based Control Rule")
#' @param xlab \code{character} X-axis label (default: "Survey Index")
#' @param ylab \code{character} Y-axis label (default: "TAC")
#' 
#' @return \code{ggplot} object
#' 
#' @examples
#' \dontrun{
#' cntrl <- FLPar(c(fmsy=0.16, btrig=5000, bbuf=2500))
#' index <- FLQuant(seq(1,10000,length.out=52), dimnames=list(year=1990:2041))
#' catch <- index %*% cntrl["fmsy"]
#' 
#' for (i in 1991:2040) {
#'   catch[,ac(i)] <- hcrSurvey(i, index, catch, cntrl, lag=0)
#' }
#' 
#' plotSurveyHCR(cntrl, index, catch, 
#'               title = "Survey-Based Control Rule",
#'               xlab = "Survey Index", 
#'               ylab = "TAC")
#' }
#' 
#' @export
plotSurveyHCR <- function(cntrl, index, catch, 
                          title = "Survey-Based Control Rule",
                          xlab = "Survey Index", 
                          ylab = "TAC") {
  
  # This function would create the visualization shown in the example
  # Implementation would depend on available plotting functions
  warning("plotSurveyHCR function not yet implemented")
  
  # Placeholder for future implementation
  # Would create a plot similar to the example in the documentation
  # showing the three zones (green, yellow, red) and the catch response
}

#' Nephrops HCR Visualization
#' 
#' @title Nephrops HCR Visualization
#' 
#' @description Creates a diagnostic plot showing the Nephrops HCR zones and catch response
#' 
#' @param cntrl \code{FLPar} Control parameters
#' @param index \code{FLQuant} Index values for plotting
#' @param catch \code{FLQuant} Catch values for plotting
#' 
#' @return \code{ggplot} object
#' 
#' @examples
#' \dontrun{
#' cntrl <- FLPar(c(fmsy=0.16, btrig=5000, bbuf=2500))
#' index <- FLQuant(seq(1,10000,length.out=52), dimnames=list(year=1990:2041))
#' catch <- index %*% cntrl["fmsy"]
#' 
#' for (i in 1991:2040) {
#'   catch[,ac(i)] <- hcrNephrops(i, index, catch, cntrl, lag=0)
#' }
#' 
#' plotNephropsHCR(cntrl, index, catch)
#' }
#' 
#' @export
plotNephropsHCR <- function(cntrl, index, catch) {
  # Call the generic survey-based plotting function with Nephrops-specific defaults
  plotSurveyHCR(cntrl = cntrl, index = index, catch = catch,
                title = "Nephrops Empirical Control Rule",
                xlab = "UWT Index", 
                ylab = "TAC")
} 