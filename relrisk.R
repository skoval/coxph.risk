# CONSIDERATION OF STARTIFICATION
projection.relrisk <- function(object, coxph.object){
	
	if(is.data.frame(x))
		return(coxph.relrisk.uncentered(coxph.object, object))
	else if (is.numeric(x))
		return(x)
	else if (class(object)=="coxph")
		return(coxph.relrisk.uncentered(coxph.object)) # DEFAULT RETURNS RELRISKS FOR DATASETS
	else
		stop(cat("No method for class ",class(object)))
	
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



