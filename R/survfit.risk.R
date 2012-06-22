# DETERMINE BASELINE SURVIVAL FROM SURVFIT
# ALLOW FOR SUBSETTING TO MAX.TIME WITH/WITHOUT STRICT INCLUSION OF MAX.TIME
baseline.survfit <- function(object, max.time, include = TRUE){
  
  sfit <- survfit(object)
  strata <- sfit$strata
  sfit.subset <- data.frame(
      surv = sfit$surv,
      time = sfit$time,
      n.event = sfit$n.event    
    )  
  sfit.subset <- subset(sfit.subset, n.event!=0)
  sfit.subset <- unique(sfit.subset)
  
  if (!is.null(strata)) {
    strata <- factor(rep(names(strata), strata), levels = names(strata))
    strata <- strata[sfit$n.event!=0]
    
    sfit.subset$haz <- sapply(levels(strata), function(level){
      h <- -log(sfit.subset$surv)[strata==level]
      c(h[1],diff(h))
    },simplify=TRUE,USE.NAMES=FALSE)
          
  }
  else{
    h <- -log(sfit.subset$surv) # CUMULATIVE BASELINE
    h <- c(h[1],diff(h)) # HAZARD INCREMENTS
    sfit.subset$haz <- h  
  }

  if(!missing(max.time)){ # SUBSET

    if(include)
      index <- sfit.subset$time<=max.time
    else
      index <- sfit.subset$time<max.time

    if(any(index))
        return(sfit.subset[index,]) 
    else # IF NO TIME BEFORE TAKE FIRST ROW ONLY
        return(sfit.subset[1,])
  }
  else
    return(sfit.subset)
}