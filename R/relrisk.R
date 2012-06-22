projection.relrisk <- function(object, data){
	
	if (is.numeric(object))
		return(object)
  else if (class(object)=="coxph")
      if (missing(data))
  		  return(coxph.relrisk.uncentered(object)) # RETURNS RELRISKS FOR FULL DATASET
	    else
	      return(coxph.relrisk.uncentered(object, data))
  else
		stop(cat("No method for class",class(object)))
	
}

# WHEN NEW DATA FRAME IS SUPPLIED NEED TO ADJUST IF STRATIFIED
coxph.relrisk.uncentered <- function(coxph.object, newdata){
	
	# RETURN VECTOR OF INDIVIDUAL RELATIVE RISKS FORM COXPH.OBJECT
	# FASTER IMPLEMENTATION THAN THROUGH FULL DESIGN
	center <- coxph.object$means%*%coef(coxph.object)
	
	if(!missing(newdata))
		lp <- predict(coxph.object, newdata, type="lp")
	else
		lp <- coxph.object$linear.predictor
		
exp(lp+center)	
}

