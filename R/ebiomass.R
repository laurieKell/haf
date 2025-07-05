#' Calculate Exploitable Biomass
#'
#' @description
#' Calculates the exploitable biomass from FLR objects using selectivity-weighted
#' catch weights and stock numbers. Exploitable biomass represents the biomass
#' that is vulnerable to fishing, weighted by the selectivity pattern.
#'
#' @param object An object of class FLStock or FLBRP
#'
#' @return An FLQuant object containing the exploitable biomass time series
#'
#' @details
#' Exploitable biomass is calculated by:
#' 1. Weighting catch weights by the selectivity pattern
#' 2. Normalizing by the peak selectivity (fapex)
#' 3. Multiplying by stock numbers
#' 4. Summing across ages
#'
#' The formula used is: \code{sum(catch.wt * selectivity / fapex(selectivity) * stock.n)}
#'
#' @examples
#' \dontrun{
#' data(ple4)
#' eb <- ebiomass(ple4)
#' plot(eb)
#' }
#'
#' @export
setGeneric("ebiomass", function(object) standardGeneric("ebiomass"))

#' @rdname ebiomass
#' @export
setMethod("ebiomass", signature(object="FLStock"),
          function(object) {
            # Get selectivity pattern
            sel <- harvest(object)
            
            # Weight catch weights by selectivity and normalize by peak
            wt <- catch.wt(object) %*% sel %/% fapex(sel)
            
            # Ensure minimum weight to avoid numerical issues
            eb.wt <- qmax(wt, 0.000001)
            
            # Calculate exploitable biomass by summing across ages
            apply(eb.wt %*% stock.n(object), 2:6, sum)
          })

#' @rdname ebiomass
#' @export
setMethod("ebiomass", signature(object="FLBRP"),
          function(object) {
            # Get selectivity pattern
            sel <- harvest(object)
            
            # Weight catch weights by selectivity and normalize by peak
            wt <- catch.wt(object) %*% sel %/% fapex(sel)
            
            # Ensure minimum weight to avoid numerical issues
            eb.wt <- qmax(wt, 0.000001)
            
            # Calculate exploitable biomass by summing across ages
            apply(eb.wt %*% stock.n(object), 2:6, sum)
          })
