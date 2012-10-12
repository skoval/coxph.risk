# CONCATENATE COXPH.OBJECTS
c.coxph.risk <- function(coxph,...){
  models <- list(coxph)
  if(!missing(..1))
    models <- c(models, list(...))
  models  
}

coxph.risk.initialize <- function(coxph,...){
  
  models <- c.coxph.risk(coxph,...)
  relrisk <- lapply(models, projection.relrisk)
  baseline <- lapply(models, survfit.risk)

  list(
    models = models,
    relrisk = relrisk,
    baseline = baseline
    )
}