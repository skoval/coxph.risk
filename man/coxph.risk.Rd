\name{coxph.risk}
\alias{coxph.risk}


\title{Absolute risk estimation with Cox proportional hazards models}

\description{
Method to estimate absolute risk in the presence of multiple competing events and with each event hazard specified by a Cox proportional hazards model.
}


\usage{
coxph.risk(interval, newdata, coxph1, ..., na.action = na.exclude)
}


\arguments{
  \item{begin}{vector specifying the begin of the projection interval, [begin, end). If scalar is given, interval is repeated for all \code{newdata}}
  \item{end}{vector specifying the end of the projection interval. If scalar is given, interval is repeated for all \code{newdata}}
  \item{newdata}{data frame containing the risk profiles for the individualized prediction as in \code{predict.glm}}
  \item{coxph1}{a \code{\link{coxph.object}} for the primary event}
  \item{...}{additional \code{coxph.objects} for competing events}
  \item{na.action}{function for handling missing data among model variables}
}

\value{
A vector of the absolute risk of the primary event occurring within [begin, end) for each individual of \code{newdata}.
}

\author{Stephanie Kovalchik <kovalchiksa@mail.nih.gov>}

%\seealso{}

\examples{

data(mgus)

# Mayo Clinic 20-35 year follow-up of patients with 
# monoclonal gammopathy of undetermined significance (MGUS)

# Hazard models of multiple myeloma, death, other plasma malignancy
# Time scale is days from MGUS diagnosis

myeloma.model <- Surv(time, status)~age+factor(sex)+alb+hgb+mspike
competing.model <- Surv(time, status)~age*factor(sex)

cox1 <- coxph(myeloma.model,data=mgus2,subset=event=="myeloma")
cox2 <- coxph(competing.model,data=mgus2,subset=event=="death")
cox3 <- coxph(competing.model,data=mgus2,subset=event=="other")

# Absolute risk predictions for multiple myeloma in 5 years
predict.data <-  mgus2[mgus2$event=="death",]

# ONLY COMPLETE CASES
predict.data <- predict.data[complete.cases(predict.data),]
risk <- coxph.risk(0, 5*365.25, newdata = predict.data,
     			     cox1, cox2, cox3)

summary(risk)

# RISK BY AGE AND GENDER AT MGUS DIAGNOSIS
cols <- c("dodgerblue","darkorchid")

plot(risk*1000~age, data = predict.data, 
	       ylab = "multiple myeloma 5-yr absolute risk (per 1000)",
	       las = 1, col = cols[predict.data$sex])

legend("topright", bty="n", levels(predict.data$sex), col= cols, pch=1)

}

\keyword{predict}

\keyword{survival}
