# BASELINE SURVIVAL/HAZARD, ALLOWING FOR STRATIFICATION
survfit.risk <- function(object){
  
  sfit <- survfit(object)
  strata <- sfit$strata
  
  if (!is.null(strata)) {
    strata <- factor(rep(names(strata), strata), levels = names(strata))

    h <- sapply(levels(strata), function(level){
      h <- -log(sfit$surv)[strata==level]
      c(h[1],diff(h))
    },simplify=TRUE,USE.NAMES=FALSE)
    
    return(data.frame(
      surv = sfit$surv,
      time = sfit$time,
      haz = as.vector(h),
      n.event = sfit$n.event,
      strata = strata
    ))         
  }
  else{
    h <- -log(sfit$surv) # CUMULATIVE BASELINE
    h <- c(h[1],diff(h)) # HAZARD INCREMENTS
  
    return(data.frame(
      surv = sfit$surv,
      time = sfit$time,
      n.event = sfit$n.event,
      haz = h
    ))
  
  }
}