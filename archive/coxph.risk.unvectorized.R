coxph.risk <- function(begin, end, projection.data, coxph1, ...){
  
  in.interval <- function(x, begin, end) 
    x$time >= begin & x$time < end
  
  models <- c.coxph.risk(coxph1, ...)
  relrisk <- sapply(models, projection.relrisk, data = projection.data)

  main <- baseline.survfit(models[[1]], max.time = end, include = FALSE)
  competing <- lapply(models[-1], baseline.survfit, max.time = main$time[nrow(main)])
  main <- main[in.interval(main, begin, end),]

  # INDIVIDUALIZED SURVIVAL FOR EACH COMPETING EVENT
  S.competing <- mapply(competing.surv,baseline=competing,rr=relrisk[-1],
                               MoreArgs=list(times=main$time))
  
  # CHECK CASE WITH SINGLE COMPETING EVENT
  pi <- main$haz*relrisk[1]*(main$surv^relrisk[1])
  pi <- cbind(pi,S.competing)
  num <- sum(apply(pi,1,prod))
  denom <- prod(c(main$surv[1]^relrisk[1],S.competing[1,]))
  
num/denom
}