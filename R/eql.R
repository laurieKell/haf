#' @title Equilibrium Model Fitting
#' @description Fits stock-recruitment models and computes biological reference points using FLR framework.
#'
#' @param object An FLStock object
#' @param model Stock-recruitment model type (default="bevholtSV")
#' 
#' @return An FLBRP object with additional attributes:
#' \itemize{
#'   \item sr - Fitted stock-recruitment model
#'   \item logLik - Model likelihood
#'   \item prod - Production characteristics
#'   \item tseries - Time series metrics
#' }
#' 
#' @examples
#' \dontrun{
#' data(ple4)
#' eql <- eql(ple4, model="rickerSV")
#' summary(eql)
#' }
#' 
#' @rdname eql
#' @export
setGeneric("eql", function(object, model) standardGeneric("eql"))

#' @rdname eql
#' @export
setMethod("eql", signature(object="FLStock"),
          eqFn<-function(object,model="bevholtSV"){
            
            spFn<-function(x){
              rfs=FLPar(c(ssb.obs(x)),dimnames=list(refpts="ssb",
                                                    quant =dimnames(refpts(x))$quant,
                                                    iter  =seq(dim(ssb.obs(x))[2])))
              rfs[,-4]=NA
              refpts(x)=rfs
              
              rtn=data.frame(model.frame(FLQuants(x,ssb=ssb.obs,catch=catch.obs),drop=TRUE),
                             sp=c(computeRefpts(x)[,"yield"]))
              rtn$pe=(c(rtn$ssb[-1]-rtn$ssb[-dim(rtn)[1]]+rtn$catch[-dim(rtn)[1]]-
                          rtn$sp[-dim(rtn)[1]],NA))/rtn$ssb
              
              rtn}
            
            spr0=spr0Yr(object)
            sr  =as.FLSR(object,model=model)
            sr  =FLCandy:::ftmb(sr,s.est =T,
                                s        =0.7, #fishlife(object)["s"],
                                s.logitsd=0.4, #fishlife(object)["sd.logit.s"],
                                spr0     =spr0)
            
            rtn=brp(FLBRP(object,nyears=dim(object)[2],
                          sr=list(model =do.call(gsub("SV","", model),list())$model,
                                  params=FLPar(apply(params(sr),1,median)))))
            
            attributes(rtn)[["sr"]]     =sr
            attributes(rtn)[["logLik"]] =logLik(sr)
            attributes(rtn)[["prod"]]   =spFn(rtn)
            attributes(rtn)[["tseries"]]=FLCandy:::tseries(object)
            attributes(rtn)[["eb.obs"]] =ebiomass(object)
            attributes(rtn)[["priors"]] =FLCandy:::tryIt(FLCandy:::calcPriors(rtn))
            attributes(rtn)[["prior2"]] =FLCandy:::tryIt(FLCandy:::getPriors(rtn))
            
            return(rtn)})

#' Simple spr0Yr function
#' @param object FLStock object
#' @return Numeric value
spr0Yr <- function(object) {
  FLCandy:::spr0Yr(object)
  }
