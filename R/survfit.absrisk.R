# NON-PARAMETRIC ABSOLUTE RISK ESTIMATION

c.survfit <- function(survfit,...){ # COLLECT SURVFIT OBJECTS
  objects <- list(survfit)
  if(!missing(..1))
    objects <- c(objects, list(...))
  objects
}

apply.survfit.risk <- function(begin, end, models){
	
	# DOES NOT EXPLICITYLY HANDLE STRATA
  
  in.interval <- function(time, begin, end) #[begin, end]
    time >= begin & time <= end
  
  event.times <- models[[1]]$time[models[[1]]$n.event>0]
  which.event.times.in.interval <- in.interval(event.times, begin, end)
  event.times.in.interval <- c(begin, event.times[which.event.times.in.interval])
    
  if(all(!event.times.in.interval))
    stop("No event times are within the projection interval.")
   
  # GET SURVIVAL TIMES FOR EACH FIT
  S <- sapply(models, function(survfit.obj){
    summary(survfit.obj, time=event.times.in.interval, extend=TRUE)$surv
  })

  # CUMULATIVE HAZARD JUST BEFORE START OF PROJECTION INTERVAL
  H0 <- -log(summary(models[[1]],time=begin*.999, extend=TRUE)$surv)
  
  nc <- length(models)
  nr <- length(event.times.in.interval)
  S <- matrix(S,nrow=nr,ncol=nc) # MAKE SURE THIS IS A MATRIX
  
  # COMPUTE HAZARD JUMPS FOR PRIMARY EVENT
  h <- -log(S[,1]) # CUMULATIVE BASELINE
  h <- c(h[1]-H0,diff(h)) # INCREMENTS
  
  # INTEGRATE
  for(i in 1:ncol(S))
    h = h*S[,i]

return(sum(h)/prod(S[1,]))
}



survfit.risk <- function(begin, end, survfit1,...){

 	fits <- c.survfit(survfit1,...)  
 
  	mapply(apply.survfit.risk, begin=begin, end=end, MoreArgs=list(models=fits))

}
  
 