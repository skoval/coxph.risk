Example of absolute risk estimation with MGUS dataset
--------------------------------------------------

Monoclonal gammopathy of undetermined significance (MGUS) is an asymptomatic condition that is 
characterized by a shift in the distribution of isotypes of immunogloblin, such that
one isotype is predominately manufactured. MGUS is believed to be an intermediate stage before
development of multiple myeloma, although most individuals diagnosed with MGUS will never develop
this cancer. Dr. Robert Kyle's study of a Mayo Clinic cohort was the first pioneering study of the
natural history of MGUS. 

In the following, we examine the outcomes for 241 of the Mayo Clinic patients, consisting of individuals
diagnosed with MGUS between the years 1956 and 1973 whose outcomes were determined 20 to 35
years after diagnosis. Our goal is to estimate the absolute risk of multiple myeloma following MGUS, 
given competing risks of death or other plasma cell malignancy. The first step of the absolute risk
calculation is to specify a relative risk model for each event type. Here, we consider Cox models and
determine risk factors among the following baseline variables as measured at the time of diagnosis: age,
gender, albumin, hemoglobin, and magnitude of monoclonal protein spike. Because a number of subjects
did not have a measured creatinine level, this measure was not considered in this set of analyses.



```r
library(survival)
```

```
## Loading required package: splines
```

```r
library(coxph.risk)

data(mgus)  # mgus2 has competing risk structure

# YEARS TO EVENT TYPE

S1 <- survfit(Surv(time, status) ~ 1, data = mgus2, subset = event == 
    "myeloma")
S2 <- survfit(Surv(time, status) ~ 1, data = mgus2, subset = event == 
    "death")
S3 <- survfit(Surv(time, status) ~ 1, data = mgus2, subset = event == 
    "other")
```




There were `39` cases of multiple myeloma, `184` deaths, and
`20` cases of other plasma malignancy during the 20 to 35 year follow-up. Because of the
risk of death, we expect that the absolute risk (crude risk or cumulative incidence)
of myeloma for intervals spanning multiple years after MGUS diagnosis could differ substantially from net risk. 

We can take a simple look at this with a small example. The net risk for myeloma within the first ten years
of diagnosis after MGUS is the estimated absolute risk if all competing events were eliminated. This can
be estimated non-parametrically by taking the product of the hazard and survival estimates over the
projection interval.



```r
end <- 365.25 * 10

# EVENT TIMES IN INTERVAL
sfit <- summary(S1)
sfit$h <- c(-log(sfit$surv[1]), diff(-log(sfit$surv)))  # hazard estimates
i <- which(sfit$time < end)

net.risk <- sum(sfit$h[i] * sfit$surv[i])
```



The net ten-year risk of multiple myeloma is `10.8`. If we now, determine the absolute risk
taking into account the competing of death and, for now, ignoring the impact of other plasma cell malignancies,
we need to include the probability of overall survival at each event time.




```r
dfit <- summary(S2, time = sfit$time[i])
crude.risk <- sum(sfit$h[i] * sfit$surv[i] * dfit$surv)
```




With a fully nonparametric analysis, results in an absolute risk of 
`8.3`. Thus the impact of death resulted in a `2.6` absolute
reduction in the risk of myeloma.

The nonparametric absolute risk can be determined directly with the `survfit.risk` function. After obtaining the primary event and competing event `survfit` objects, the absolute risk estimate for a given time interval [begin, end) can be obtained as follows:



```r
survfit.risk(0, 5 * 365.25, S1, S2, S3)
```

```
## [1] 0.03728
```




To specify multiple event times, use vector arguments for the `begin` and `end`. In what follows, we look at the absolute risk of death within 1 to 10 years of MGUS diagnosis in yearly intervals.



```r
survfit.risk(rep(0, 10), 1:10 * 365.25, S2, S1, S3)
```

```
##  [1] 0.05424 0.08760 0.12080 0.14124 0.18937 0.21311 0.25541 0.28152
##  [9] 0.33213 0.35688
```




For semiparametric estimates of risk, where we can introduce the impact of risk factors through Cox's proportional hazards model, we would use the function `coxph.risk`. The syntax is similar to `survfit.risk` only the modeling objects are of the `coxph` class.



```r
base.model <- Surv(time, status) ~ age + factor(sex) + alb + hgb + 
    mspike

cox1 <- coxph(base.model, data = mgus2, subset = event == "myeloma")
cox2 <- coxph(base.model, data = mgus2, subset = event == "death")
cox3 <- coxph(base.model, data = mgus2, subset = event == "other")

coxph.risk(0, end, newdata = mgus2[1, ], cox1, cox2, cox3)
```

```
## [1] 0.06694
```



