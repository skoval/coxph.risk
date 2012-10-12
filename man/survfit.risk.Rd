\name{survfit.risk}
\alias{survfit.risk}


\title{Non-parametric absolute risk estimation}

\description{
Method to estimate a fully non-parametric absolute risk in the presence of multiple competing events. 
}


\usage{
survfit.risk(begin, end, survfit1, ...)
}


\arguments{
  \item{begin}{vector with each row specifying the beginning time in the [begin, end) of the projection interval}
  \item{end}{vector with each row specifying the ending time in the [begin, end) of the projection interval}
  \item{survfit1}{a \code{\link{survfit.object}} for the primary event}
  \item{...}{additional \code{survfit.objects} for competing events}
}

\value{
A vector of the absolute risk of the primary event occurring within [begin, end).
}

\details{
The non-parametric approch assumes that the time-to-events of all individuals in the sample in which the event times were observed follow the same distribution. This is equivalent to assuming that the hazard ratio for any two individuals of the sample is 1 for all time points.
}

\author{Stephanie Kovalchik <kovalchiksa@mail.nih.gov>}

\seealso{\code{\link{coxph.risk}}}

\examples{

data(mgus)

# Mayo Clinic 20-35 year follow-up of patients with 
# monoclonal gammopathy of undetermined significance (MGUS)

# Hazard models of multiple myeloma, death, other plasma malignancy
# Time scale is days from MGUS diagnosis

surv1 <- survfit(Surv(time, status)~1,data=mgus2,subset=event=="myeloma")
surv2 <- survfit(Surv(time, status)~1,data=mgus2,subset=event=="death")
surv3 <- survfit(Surv(time, status)~1,data=mgus2,subset=event=="other")

# Absolute risk predictions for multiple myeloma consecutive biyearly intervals
bytwo <- seq(0,9,by=2)
intervals <- cbind(bytwo, bytwo+2)*365.25
risk <- survfit.risk(intervals[,1],intervals[,2], surv1, surv2, surv3)
risk

# How much do competing risks of death and other plasma malignancy reduce net risk?
netrisk <- survfit.risk(intervals[,1], intervals[,2], surv1)
netrisk

names(risk) <- paste("[",bytwo,", ",bytwo+2,")",sep="")

dotchart(risk*100, xlim = range(c(risk, netrisk)*100), 
                    xlab = "risk of multiple myeloma progression (\%)",
                    ylab = "projection interval")
                    
points(x=netrisk*100, y=1:length(risk), col="red")

legend("topleft",bty="n",c("Net Risk", "Abs. Risk"), pch=1,col = c("red", "black"))


# How does this compare with cumulative risk?
intervals <- cbind(rep(0, length(bytwo)), bytwo+2)*365.25
risk <- survfit.risk(intervals[,1], intervals[,2], surv1, surv2, surv3)
netrisk <- survfit.risk(intervals[,1], intervals[,2], surv1)

names(risk) <- paste("[0, ",bytwo+2,")",sep="")

dotchart(risk*100, xlim = range(c(risk, netrisk)*100), 
                    xlab = "risk of multiple myeloma progression (\%)",
                    ylab = "projection interval")
                    
points(x=netrisk*100, y=1:length(risk), col="red")

legend("topleft",bty="n",c("Net Risk", "Abs. Risk"), pch=1,col = c("red", "black"))

### STRATIFIED BY GENDER

surv1 <- survfit(Surv(time, status)~1,data=mgus2,subset=event=="myeloma"&sex=="male")
surv2 <- survfit(Surv(time, status)~1,data=mgus2,subset=event=="death"&sex=="male")
surv3 <- survfit(Surv(time, status)~1,data=mgus2,subset=event=="other"&sex=="male")

# 5-year, 10-year RISK OF MULTIPLE MYELOMA
intervals <- cbind(c(0,0), c(5,10))*365.25
male.risk <- survfit.risk(intervals[,1], intervals[,2], surv1, surv2, surv3)

male.risk

}

\keyword{predict}

\keyword{survival}
