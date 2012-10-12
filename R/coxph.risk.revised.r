basehaz.coxph.risk.revised <- function(object, times){
  
  sfit <- summary(survfit(object), time=times, extend=TRUE)
  H <- -log(sfit$surv)
  z0 <- object$means
  bz0 <- sum(z0 * coef(object))
  H <- H * exp(-bz0)
  S <- exp(-H) # BASELINE
 
  # WHEN START OF PROJECTION IS NOT BEFORE FIRST EVENT
  # 	ADJUST FIRST HAZARD JUMP; ONLY IMPACTS PRIMARY EVENT
  
  H0 <- -log(summary(survfit(object), time=times[1]*.999, extend=TRUE)$surv)
  h <- c(H[1]-H0,diff(H))
  
data.frame(time=times,haz=h,surv=S/S[1])
}

coxph.risk.baseline <- function(begin,end,models){

 in.interval <- function(x, begin, end) 
    x >= begin & x <= end
  
  # TAKE THE SET OF EVENT TIMES AND RETURN THE CUMULATIVE HAZARDS FOR EACH MODEL
  # IF NO TIMES, THE HAZARD IS ZERO
  event.times <- survfit(models[[1]])
  event.times <- event.times$time[event.times$n.event>=1]
  which.event.times <- in.interval(event.times, begin, end)
  
    if(all(!which.event.times))
	  	NA
	else{

	  event.times <- event.times[which.event.times]

	  Hazards <- mapply(basehaz.coxph.risk.revised, object=models,
	  						MoreArgs=list(times=event.times),SIMPLIFY=FALSE)
  	  Hazards
  }
}



coxph.risk <- function(begin, end, newdata, coxph1, ...){
 
  risk <- function(begin, end, models, RR){ 
  	
  	H <- coxph.risk.baseline(begin, end, models) # CALLS BASELINE ONCE
   	if(!is.list(H))
  		0
  	else{
  		absrisk <- H[[1]]$surv^RR[1]*H[[1]]$haz*RR[1]
  		for(i in 2:length(H)){
  			absrisk <- absrisk*(H[[i]]$surv^RR[i])
  		}
  	sum(absrisk)
  	}
  }
  
  risk.fixed.interval <- function(H,RR){ 
   	if(!is.list(H))
  		0
  	else{
  		absrisk <- H[[1]]$surv^RR[1]*H[[1]]$haz*RR[1]
  		for(i in 2:length(H)){
  			absrisk <- absrisk*(H[[i]]$surv^RR[i])
  		}
  	sum(absrisk)
  	}
  }
  
  models <- c.coxph.risk(coxph1, ...) # COLLECT COX MODELS
  AllVars <- unique(unlist(sapply(models, function(x) all.vars(x$formula))))
  newdata <- subset(newdata, select = unique(AllVars))
  which.kept <- complete.cases(newdata)
  
  if(!all(which.kept)){
  	warning("Missing cases excluded.")
  	
  	if(length(begin)>1){
	  	begin <- begin[which.kept]
  		end <- end[which.kept]
  	}
  	
  	newdata <- newdata[which.kept,]
  }
  
  # CAUSE-SPECIFIC RR, VECTOR FOR MULTIPLE PROJECTIONS
  # NOT LIKING THIS IMPLEMENTATION
  rr <- sapply(models, projection.relrisk, data = newdata) # MULTIPLE ROWS FOR EACH EVENT TYPE

  if(is.matrix(rr))
	  rr.list <- lapply(1:nrow(rr), function(x) rr[x,]) # RETURN LIST BY EACH ROW
  else
  	  rr.list <- list(rr)
	
  if(length(begin)==1){
	  H <- coxph.risk.baseline(begin, end, models) # CALLS BASELINE ONCE  	
	  risks <- mapply(risk.fixed.interval,RR=rr.list,MoreArgs=list(H=H))
  }
  else{
	  risks <- mapply(risk,begin=begin,end=end,RR=rr.list,MoreArgs=list(models=models))
}

risks
}
	 
 