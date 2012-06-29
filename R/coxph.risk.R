coxph.risk <- function(interval, newdata, coxph1, ..., na.action = na.exclude){
	  
  begin <- interval[1]
  end <- interval[2]
  
  in.interval <- function(x, begin, end) 
    x$time >= begin & x$time < end
  
  models <- c.coxph.risk(coxph1, ...)
  AllVars <- unique(unlist(sapply(models, function(x) all.vars(x$formula))))
  newdata <- subset(newdata, select = unique(AllVars))
  newdata <- na.action(newdata)
  
  # REMOVING MISSING VARIABLES IN NEW DATA
  relrisk <- lapply(models, projection.relrisk, data = newdata) # MULTIPLE ROWS FOR EACH EVENT TYPE
  
  main <- baseline.survfit(models[[1]], max.time = end, include = FALSE)
  competing <- lapply(models[-1], baseline.survfit, max.time = main$time[nrow(main)])
  main <- main[in.interval(main, begin, end),]

  # CREATE MATRICES FOR MULTIPLE RISK COMPUTATION
  # EACH COLUMN IS A UNIQUE RISK TYPE
  nc = nrow(newdata)
  nr = length(main$haz)  
  relrisk <- lapply(relrisk, function(x) matrix(x, nr, nc, byrow=TRUE))
  
  # MATRIX FORM, REPEAT HAZARD BY ROW: REPEAT OPERATIONS BY ROW
  H <- matrix(main$haz, nr, nc)
  S0 <- matrix(main$surv, nr, nc)^relrisk[[1]]
  haz0 <- H*relrisk[[1]]
  
  # INDIVIDUALIZED SURVIVAL FOR EACH COMPETING EVENT
  S.competing <- mapply(competing.surv,baseline=competing,SIMPLIFY=FALSE,
  											MoreArgs=list(times=main$time, rr=1))
 
  S.competing <- lapply(S.competing, function(x) matrix(x, nr, nc))
    
  for(i in 1:length(S.competing)){
  	S0 <- S0*S.competing[[i]]^relrisk[[i+1]]
  }
  
  denom <- S0[1,] # SURVIVAL TO START 
  num <- colSums(S0*haz0)
  risk <- num/denom
  names(risk) <- row.names(newdata)
  
 risk
}