#' Fits Stock Recruitment Relationships (SRR) in TBM (with flexible priors for s and r0)
#'
#' @param object Input FLSR object.
#' @param s steepness parameter of SRR (fixed or prior mean)
#' @param spr0 unfished spawning biomass per recruit from FLCore::spr0(FLStock) or FLCore:::spr0Yr
#' @param s.est option to estimate steepness
#' @param s.logitsd prior sd for logit(s), default is 1.3 (flat) if s.est = TRUE
#' @param inflect Inflection point for the segreg model. If NA (default), it is estimated; if numeric, it is fixed at the provided value.
#' @param inits option to specify initial values of log(r0), log(SigR) and logit(s)
#' @param lower option to specify lower bounds of log(r0), log(SigR) and logit(s)
#' @param upper option to specify upper bounds of log(r0), log(SigR) and logit(s)
#' @param SDreport option to converge hessian and get vcov
#' @param prior_s Optional prior mean for steepness (s)
#' @param cv_s Optional coefficient of variation for steepness prior
#' @param prior_r0 Optional prior mean for r0
#' @param cv_r0 Optional coefficient of variation for r0 prior
#'
#' @return A list containing elements 'FLSR', of class *FLSR*
#'
#' @examples
#' # See ftmb for usage
ftmb2<-function(object, 
               spr0=spr0, 
               model=SRModelName(FLCore:::model(object)),
               s=0.7, s.est=TRUE, s.logitsd=1.3,
               inflect=NA,
               inits=function(object=object,s=s,inflect=inflect) {
                 inflect_init <- ifelse(is.na(inflect), 0, log(inflect))
                 c(median(c(log(rec(object)))), log(0.40), to_logits(s), inflect_init)
               },
               lower=function(object=object,s=-20,inflect=inflect) {
                 inflect_lower <- ifelse(is.na(inflect), -20, log(inflect))
                 c(quantile(c(log(rec(object))),probs=0.1), log(0.05), -20, inflect_lower)
               },
               upper=function(object=object,s= 20,inflect=inflect) {
                 inflect_upper <- ifelse(is.na(inflect), 20, log(inflect))
                 c(quantile(c(log(rec(object))),probs=0.9), log(1.50), 20, inflect_upper)
               },
               short   =TRUE,
               SDreport=short,
               prior_s=NULL, cv_s=NULL, prior_r0=NULL, cv_r0=NULL,
               return_iters=FALSE) {
  
     # Check if data has iterations
   rec_dims <- dim(rec(object))
   ssb_dims <- dim(ssb(object))
   has_iters <- length(rec_dims) > 2 && rec_dims[6] > 1
   
       # Check for dimension consistency between rec and ssb
    # This can happen when min age is not 0 - rec has fewer year classes than ssb
    if(has_iters && !identical(rec_dims, ssb_dims)) {
      # Check if this is the expected year class mismatch
      if(rec_dims[2] < ssb_dims[2] && rec_dims[6] == ssb_dims[6]) {
        # This is the expected case where rec has fewer year classes than ssb
        # We can still proceed with iterations, but need to handle the dimension mismatch
        warning("Recruitment has fewer year classes than SSB (likely due to age structure). Proceeding with iteration handling.")
      } else {
        # Unexpected dimension mismatch
        warning("SSB and recruitment have different dimensions. This may be due to age structure issues. Proceeding with single iteration fit.")
        has_iters <- FALSE
      }
    }
  
  # If no iterations, proceed with original logic
  if(!has_iters) {
    # Original ftmb2 logic for single iteration
    if(is.null(s)& s.est) s=0.6 # central value
    #if(is.null(s)&!s.est){s=0.8}
    
    # GET rec, ssb
    rec=c(rec(object))
    ssb=c(ssb(object))
    
    if ("FLQuant"%in%is(spr0))
      spr0.=c(spr0)
    else
      spr0.=rep(spr0,length(rec))
    
    # SET init and bounds
    inits.=inits(object,s,inflect)
    lower.=lower(object,inflect=inflect)
    upper.=upper(object,inflect=inflect)
    
    # --- PRIORS ---
    # Steepness prior: default tologit(s), s.logitsd unless prior_s/cv_s provided
    if(!is.null(prior_s) & !is.null(cv_s)) {
      prior_s_logit <- to_logits(prior_s)
      prior_s_sd <- sqrt(log(1 + cv_s^2))
      prior_s_vec <- c(prior_s_logit, prior_s_sd)
    } else {
      prior_s_vec <- c(to_logits(s), s.logitsd)
    }
    # r0 prior: default NULL unless prior_r0/cv_r0 provided
    if(!is.null(prior_r0) & !is.null(cv_r0)) {
      prior_r0_log <- log(prior_r0)
      prior_r0_sd <- sqrt(log(1 + cv_r0^2))
      prior_r0_vec <- c(prior_r0_log, prior_r0_sd)
    } else {
      prior_r0_vec <- NULL
    }
    
    # SET TMB input
    inp=list(
      # data
      Data = list(ssb=ssb, rec=rec, prior_s = prior_s_vec, spr0=spr0., nyears=length(ssb),
                  # model
                  Rmodel = which(model==c("bevholtSV","rickerSV","segreg"))-1),
                  # inits
                  Params = list(log_r0=inits.[1], log_sigR=inits.[2],logit_s=inits.[3], log_inflect=inits.[4]),
                  # bounds
                  lower=lower., upper=upper.,
                  #
                  ReportSD = SDreport)
    # Always add prior_r0 (empty if not used)
    inp$Data$prior_r0 <- if(!is.null(prior_r0_vec)) prior_r0_vec else numeric(0)
    
    # Compile TMB inputs 
    Map=list()
    # Turn off steepness estimation
    if(!s.est) Map[["logit_s"]] = factor( NA ) 
    # Turn off inflection estimation if inflect is not NA
    if(!is.na(inflect)) Map[["log_inflect"]] = factor( NA )
    
    # CREATE TMB object
    Obj=TMB::MakeADFun(data=inp$Data, parameters=inp$Params, map=Map, DLL="FLCandy", silent=TRUE)
    
    capture.output({Opt=stats::nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
                                      control=list("trace"=1, "eval.max"=1e4, "iter.max"=1e4),
                                      lower=inp$lower, upper=inp$upper)})
    
    Opt[["diagnostics"]] = data.frame(Est=Opt$par,final_gradient=Obj$gr(Opt$par))
    
    Report=Obj$report()
    
    if(SDreport) SD=try(TMB::sdreport(Obj))
    
    # LOAD output in FLSR
    
    # DEBUG HACK
    model(object)=switch(model, bevholtSV=bevholt, rickerSV=ricker,segreg=segreg)
    
    fitted(object)=c(Report$rec_hat)
    residuals(object)=log(rec(object)) - log(fitted(object))
    
    #params(object)=FLPar(a=Report$a, b=Report$b)
    
    params(object)=FLPar(s=Report$s, R0=Report$r0)
    params(object)=propagate(params(object),length(spr0.))
    params(object)=rbind(params(object),FLPar(spr0=spr0.))
    params(object)=rbind(params(object),FLPar(v=params(object)["R0"]*params(object)["spr0"]))[c("s","v","spr0")]
    
    # For segreg model, use the inflection point directly as parameter b
    if(model=="segreg") {
      # Get the inflection point from the TMB fit or use the fixed value
      if(is.na(inflect)) {
        # Inflection point was estimated
        inflect_val = exp(Report$log_inflect)
      } else {
        # Inflection point was fixed
        inflect_val = inflect
      }
      # For segreg: a = r0/b (the slope), b = inflection point
      # The formula is: rec = a * min(ssb, b)
      a_val = Report$r0 / inflect_val
      params(object)=FLPar(a=a_val, b=inflect_val)
    } else {
      # For other models, use the standard ab() conversion
      par=as(as.data.frame(aaply(params(object)[c("s","v","spr0")],2,
                                 function(x) ab(FLPar(x),gsub("SV","",model)))),"FLPar")[c("a","b")]
      params(object)=par
    }
    
    # Always return the FLSR object, regardless of short parameter
    # Ensure we have all required slots
    if(!is.null(Opt)) {
      logLik(object)=-Opt$objective
      logLik(object)["df"]=df(s.est, inflect)
    }
    return(object)
      
    df<-function(s.est, inflect) {
      # Base number of parameters
      num_params <- 2  # log_r0 and log_sigR are always estimated
      # Add logit_s if it is being estimated
      if (s.est) num_params <- num_params + 1
      # Add log_inflect if it is being estimated
      if (is.na(inflect)) num_params <- num_params + 1
      return(num_params)}
    
    logLik(object)=-Opt$objective
    logLik(object)["df"]=df(s.est, inflect)
    
    return(object)
  }
  
     # Handle iterations
   else {
     # Get dimensions
     n_iters <- rec_dims[6]
     
     # Fit model to each iteration
     results <- lapply(1:n_iters, function(i) {
       
       # Create subset for this iteration
       object_iter <- object
       rec(object_iter) <- rec(object)[,,,,,i]
       ssb(object_iter) <- ssb(object)[,,,,,i]
       
       # Handle spr0 for this iteration
       spr0_iter <- if("FLQuant" %in% is(spr0)) spr0[,,,,,i] else spr0
       
       # Fit model and extract parameters
       tryCatch({
         result <- ftmb2(object_iter, spr0=spr0_iter, SDreport=FALSE, 
                        model=model, s=s, s.est=s.est, s.logitsd=s.logitsd, inflect=inflect,
                        prior_s=prior_s, cv_s=cv_s, prior_r0=prior_r0, cv_r0=cv_r0)
         # Extract parameters from FLSR object
         c(params(result))
       }, error = function(e) {
         NULL
       })
     })
    
    # Check if any fits succeeded
    successful_fits <- !sapply(results, is.null)
    if(!any(successful_fits)) {
      stop("No successful fits across iterations")
    }
    
         # Get parameter names from first successful fit
     param_names <- names(results[[which(successful_fits)[1]]])
     n_params <- length(param_names)
     
     # Create parameter array - just store the parameter values, not year dimension
     param_array <- array(NA, dim=c(n_params, n_iters),
                         dimnames=list(params=param_names, 
                                      iter=1:n_iters))
     
     # Fill array with results
     for(i in which(successful_fits)) {
       param_values <- results[[i]]
       param_array[,i] <- param_values
     }
    
         # Create an FLSR object with the fitted parameters
     # Use median parameters for the main object
     median_params <- apply(param_array, 1, median, na.rm=TRUE)
     
     # Create FLSR object with median parameters
     result_object <- object
     params(result_object) <- FLPar(median_params)
     
           # For failed iterations, we need to ensure the FLSR object has consistent dimensions
      # Use only the successful iterations for the data slots
      successful_iter_indices <- which(successful_fits)
      
      # Update rec and ssb to only include successful iterations
      if(length(successful_iter_indices) < n_iters) {
        rec(result_object) <- rec(object)[,,,,,successful_iter_indices, drop=FALSE]
        ssb(result_object) <- ssb(object)[,,,,,successful_iter_indices, drop=FALSE]
      }
      
      # Ensure the FLSR object has consistent dimensions
      # The fitted and residuals should match the recruitment dimensions (fewer year classes)
      rec_dims_final <- dim(rec(result_object))
      ssb_dims_final <- dim(ssb(result_object))
      
      # If there's still a dimension mismatch, use the recruitment dimensions as the reference
      if(rec_dims_final[2] != ssb_dims_final[2]) {
        # Use the recruitment year range for the final object
        # This ensures the FLSR object has consistent dimensions
        ssb(result_object) <- ssb(result_object)[,1:rec_dims_final[2],,,,, drop=FALSE]
      }
     
     # Ensure we have fitted and residuals slots
     # Create fitted values using median parameters
     fitted_vals <- predict(result_object, ssb=ssb(result_object))
     fitted(result_object) <- fitted_vals
     residuals(result_object) <- log(rec(result_object)) - log(fitted_vals)
    
    # Set a dummy logLik (since we don't have the full optimization results)
    logLik(result_object) <- structure(-999, df=length(median_params), class="logLik")
    
    # Ensure the model is set correctly
    model(result_object) <- switch(model, bevholtSV=bevholt, rickerSV=ricker, segreg=segreg)
    
         # If return_iters is TRUE, add the full iteration results as an attribute
     if(return_iters) {
       # Convert array to FLPar with proper dimensions
       iter_params <- FLPar(param_array)
       attr(result_object, "iterations") <- iter_params
     }
    
    return(result_object)
  }
} 