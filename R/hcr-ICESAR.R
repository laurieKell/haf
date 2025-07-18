#' ICES Advice Rule (Hockey-Stick HCR)
#' 
#' @title ICES Advice Rule (Hockey-Stick HCR)
#' 
#' @description Implements the ICES advice rule (hockey-stick harvest control rule) to calculate Total Allowable Catch (TAC) based on stock biomass relative to reference points. This is the standard ICES approach for category 1 and 2 stocks.
#' 
#' @author Laurence Kell, Sea++
#'  
#' @param object An object of class \code{FLStock} 
#' @param eql \code{FLBRP} with a stock recruitment relationship used for projection
#' @param sr_deviances \code{FLQuant} recruitment deviates on the log scale (multiplicative)
#' @param params \code{FLPar} HCR parameters, specifying blim, btrig, bmin, ftar and fmin
#' @param start \code{numeric} first year for simulation (default: max year - 10)
#' @param end \code{numeric} last year for simulation (default: start + 10)
#' @param interval \code{numeric} time step, 1 year by default
#' @param lag \code{numeric} lag between assessment and advice (default: 1)
#' @param err \code{FLQuant} assessment error on SSB for year used for HCR
#' @param implErr \code{numeric} or \code{FLQuant} implementation error
#' @param bndTac \code{numeric} bounds on TAC, by default these are turned off, for 20 percent constraint set to c(0.8,1.2)
#' @param bndWhen \code{character} when to apply bounds, default "btrig"
#' @param bndCap \code{numeric} bound cap for large changes
#' @param perfect \code{logical} whether to use perfect knowledge
#' @param ... any additional arguments
#' 
#' @return Returns a \code{list} with \code{FLStock} and \code{FLPar} objects for the stock and HCR 
#'
#' @details
#' The hockey-stick HCR implements the standard ICES approach where:
#' \itemize{
#'   \item F = Ftarget when B >= Btrigger
#'   \item F = Fmin when B <= Blim
#'   \item F decreases linearly between Blim and Btrigger
#' }
#' 
#' The rule includes optional TAC bounds and implementation error to simulate real-world conditions.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ple4)
#' 
#' # Set up HCR parameters
#' params <- FLPar(
#'   blim = 100000,   # Biomass limit reference point
#'   btrig = 150000,  # Biomass trigger reference point  
#'   bmin = 100000,   # Minimum biomass
#'   ftar = 0.3,      # Target fishing mortality
#'   fmin = 0.05      # Minimum fishing mortality
#' )
#' 
#' # Run HCR simulation
#' result <- hcrICESAR(ple4, eql, sr_deviances, params, 
#'                   start = 2010, end = 2030)
#' }
#' 
#' @export
setGeneric('hcrICESAR', function(object,eql,...) standardGeneric('hcrICESAR'))

setMethod('hcrICESAR', signature(object="FLStock",eql='FLBRP'),
          function(object,eql,sr_deviances,params,
                   start   =max(dimnames(object)$year)-10,
                   end     =start+10,
                   interval=1,
                   lag     =1,
                   err     =NULL,
                   implErr =0,
                   bndTac  =c(0,Inf),
                   bndWhen ="btrig",
                   bndCap  =1e6,...){

            if ("perfect"%in%names(list(...)))
              perfect=list(...)$perfect
            else 
              perfect=FALSE

            bias<-function(object,eql,err,iYr,lag=1){
              stock.n(object)[,ac(iYr)]=stock.n(object)[,ac(iYr)]%*%err[,ac(iYr)]
              
              if (lag>0)
                object=fwd(object,catch=catch(object)[,ac(iYr-(lag:0))],sr=eql)

              window(object,end=iYr+1)}
           
            chk=NULL
            ## Loop round years
            cat('\n==')
            for (iYr in seq(start,end,interval)){
              cat(iYr,", ",sep="")
              
              stkYrs=iYr-lag
              refYrs=iYr
              hcrYrs=iYr+1
              
              #refpts(eql)[]=1
              #refpts(eql)=propagate(refpts(eql),dim(object)[6])

              #if (!is.null(err))
              #  mp=bias(object,eql,err,iYr,lag=hcrYrs-stkYrs-1)
              #else
                mp=object
  
              if ("FLPar"%in%is(params)) {
                par=params[,min(iYr-start+1,dims(params)$iter),drop=T]}
              else
                par=params
  
              res=hcrFn(mp,eql,par,
                        stkYrs,
                        refYrs,
                        hcrYrs,
                        perfect=perfect,
                        sr_deviances=sr_deviances)
              
              ##TAC Bounds, only if stock above lower ref point
              flag=c(ssb(mp)[,ac(stkYrs)]>c(FLPar(params)[bndWhen])[1])
              rtn =res[[1]]

              ### bound TAC
              rtn[,,,,,flag]=qmax(rtn,catch(object)[,ac(refYrs)]*bndTac[1])[,,,,,flag]
              rtn[,,,,,flag]=qmin(rtn,catch(object)[,ac(refYrs)]*bndTac[2])[,,,,,flag]
              
              ## if change > bndCap then no bound
              if (bndCap<1){ 
                ctcRatio=res[[1]]%/%(catch(object)[,ac(refYrs)])
                flag=flag&(ctcRatio<(1-bndCap))|(ctcRatio>(1+bndCap))
                
              if (any(flag))
                rtn[,,,,,flag]=res[[1]][,,,,,flag]
              }

              ## Implementation Error
              if ("FLQuant"%in%is(implErr))
                rtn=rtn%*%(1+implErr[,dimnames(rtn)$year])
   
              object=fwd(object,catch=rtn,sr=eql,residuals=sr_deviances)
              chk=rbind(chk,data.frame(iter=factor(seq(dim(rtn)[6]),
                                               levels=seq(seq(dim(rtn)[6])),
                                               labels=seq(dim(rtn)[6])),
                                       res[[2]],catch=c(rtn[,1])))
            }
            cat("==\n")
            
        list(window(object,end=end),chk)})

#' @title Internal HCR Function
#' @description Internal function implementing the core hockey-stick HCR logic
#' @keywords internal
hcrFn<-function(object,eql,params,
                stkYrs=max(as.numeric(dimnames(stock(object))$year))-2,
                refYrs=max(as.numeric(dimnames(catch(object))$year))-1,
                hcrYrs=max(as.numeric(dimnames(stock(object))$year)),
                maxF  =2,
                nyrs  =3,
                sr_deviances=NULL,
                perfect=FALSE,
                ...) {

  
  ## HCR
  setTarget<-function(params,stk,hcrYrs){
    
    if (length(dim(params))==2){
      flag=c(1,2)[1+as.numeric(c(stk)>c(params["btrig","lower"]))]
    }else{
      params=propagate(FLPar(params),dim(stk)[6])
      flag=rep(1,dim(stk)[6])}
    
    a=FLPar(a=array((params['ftar',flag]-params['fmin',flag])/(params['btrig',flag]-params['bmin',flag]),c(1,dim(stk)[6])))
    b=FLPar(b=array( params['ftar',flag]-a*params['btrig',flag],c(1,dim(stk)[6])))    
    rtn=(stk%*%a)  
    rtn=FLCore::sweep(rtn,2:6,b,'+')
    
    fmin=FLQuant(c(params['fmin',flag]),dimnames=list(iter=seq(dim(stk)[6])))
    ftar=FLQuant(c(params['ftar',flag]),dimnames=list(iter=seq(dim(stk)[6])))
    
    for (i in seq(dims(object)$iter)){
      FLCore::iter(rtn,i)[]=max(FLCore::iter(rtn,i),FLCore::iter(fmin,i))
      FLCore::iter(rtn,i)[]=min(FLCore::iter(rtn,i),FLCore::iter(ftar,i))} 
    
    rtn=window(rtn,end=max(hcrYrs))
    if (hcrYrs>(stkYrs+1))
      rtn[,ac((stkYrs+1):(hcrYrs-1))]=fbar(object)[,ac((stkYrs+1):(hcrYrs-1))] #####bug
    for (i in hcrYrs)
      rtn[,ac(i)]=rtn[,dimnames(rtn)$year[1]]
    
    chk=data.frame(hcrYrs=hcrYrs,
                   ssb   =c(stk),
                   f     =c(rtn[,ac(hcrYrs)]))
    list(rtn[,ac(hcrYrs)],chk)}
  
  if (!perfect)
    status=FLCore::apply(ssb(object)[,ac(stkYrs)],6,mean)
  else  
    status=FLCore::apply(ssb(object)[,ac(min(hcrYrs))],6,mean)
  
  res=setTarget(params,status,hcrYrs)
  rtn=res[[1]]
  
  hvt=rtn
  
  ## TACs for target F
  object=window(object, end=max(as.numeric(hcrYrs)))
  
  ## short term projection setting of future vectors
  if (!perfect){
    ## set up vectors for in year projection and hcr
    for (i in slotNames(FLStock())[-c(1,4,7,10,18:20)])
      slot(object,i)[,ac(min(hcrYrs))]=apply(slot(object,i)[,ac(stkYrs-(nyrs)-1)],c(1,6),mean)
  
    if (length(hcrYrs)>1)
      for (i in slotNames(FLStock())[-c(1,4,7,10,18:20)])
        for (j in hcrYrs[-1])
          slot(object,i)[,ac(j)]=slot(object,i)[,ac(hcrYrs[1])]

    if (stkYrs<(hcrYrs-1))
       object=fwd(object,fbar=fbar(object)[,ac(max(stkYrs+1):min(as.numeric(hcrYrs)-1))],sr=eql)
    }
  
  hvt[is.na(hvt)]=6.6666
  
  if (!perfect)
    rtn=catch(fwd(object, fbar=hvt,sr=eql))[,ac(hcrYrs)]
  else  
    rtn=catch(fwd(object, fbar=hvt,sr=eql,residuals=sr_deviances))[,ac(hcrYrs)]
  rtn[]=rep(c(apply(rtn,c(3:6),mean)),each=dim(rtn)[2])
  
  return(list(rtn,res[[2]]))} 