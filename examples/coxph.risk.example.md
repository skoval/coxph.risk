Absolute risk formulation
========================================================

The definition of absolute risk is the probability that a disease-free individual will develop disease within a specified interval of time before any competing event given their individual risk factors. Letting $m$ ($1,\dots,M$) index the number of failure types, the formula for the absolute risk of experiencing the $m$ th event within the time interval $[t_0, t_1)$ in the presence of $M-1$ competing events is

$$
\pi(t_0,t_1;\vec{x}) = \left[ \prod_{i=1}^M S_{i}(t_0;\vec{x}_i)  \right]^{-1} \int_{t_0}^{t_1} \lambda_{m}(u;\vec{x}_m) \prod_{i=1}^M S_{i}(u;\vec{x}_i) du
$$

where $\vec{x} = (\vec{x}_1,\dots,\vec{x}_M)$ is a set of  cause-specific covariate vectors, $S_m(u;x_m) = \exp (- \int_0^u \lambda_m(v;\vec{x}_m) dv )$ and $\lambda_{m}(u;\vec{x}_m)$ are the cause-specific survival and hazard functions given covariates $\vec{x}_m$ [(Benichou 1990)](http://www.ncbi.nlm.nih.gov/pubmed/2242416).

We suppose that each event type follows the semiparametric model

$$
\lambda_m(v;\vec{x}_m) = \lambda_{0m}(v) r(\beta_m,\vec{x}_m)
$$

with the arbitrary baseline hazard function $\lambda_{0m}(v)$ and relative risk model that is a function of the covariates $\vec{x}_m$ denoted as $r(\beta_m,\vec{x}_m)$. Both of these components are strictly positive over all time.

When no distributional assumption is made for $\lambda_{0m}$ the estimator for the cause-specific risk of the primary event ($i = 1$) within the interval $[t_0, t_1)$, given $\vec{x}$, is

$$
\hat{\pi}(t_0,t_1;\vec{x}) = \left[ \prod_{i=1}^M \hat{S}_{0i}(t_0)^{r (\hat{\beta}_i, \vec{x}_i )}  \right]^{-1}  r (\hat{\beta}_1, \vec{x}_1 )\sum_{t_0 \leq u < t_1} \hat{\lambda}_{01}(u) \prod_{i=1}^M \hat{S}_{0i}(u)^{r(\hat{\beta}_i, \vec{x}_i )},
$$

where $\hat{S}_{0i}(u)$ is the cause-specific baseline survival function and $\hat{\lambda}_{01}(u)$ the primary-event baseline hazard function at time $u$. 

When the Cox proportional hazard is used for the relative risk model, each of the components for the absolute risk can be readily computed from most standard software packages. The `coxph.risk` puts all of the steps together for this case using a  Nelson-Aalen estimator[(Aalen 1978)](http://www.jstor.org/stable/10.2307/2958557) for the cause-specific baseline hazard function of each failure type. 

***

Example: Cox-Exp model
--------------------------------------------------------

Suppose the disease of interest has one competing event. Let the time to disease be $T_1$ and the time to the competing event be $T_2$. As a test model, we will suppose that, conditional on a set of respective risk factors $x_1$ and $x_2$, $T_1$ and $T_2$ are independent exponential random variables, whose hazards follow Cox's proportional hazards model.

The distribution for $T_1$ is

$$
T_1 | x_1 \sim \mbox{Exp}(\lambda_{01} \exp (\beta'x_1) )
$$

and $T_2$ is

$$
T_2 | x_2 \sim \mbox{Exp}(\lambda_{02} \exp (\gamma'x_2) ).
$$

The absolute risk of disease between time $\tau_0$ and $\tau_1$ is

$$
\pi = \lbrack S_1 (\tau_0) S_2 (\tau_0) \rbrack^{-1} \int_{\tau_0}^{\tau_1} \lambda_{1} (u) S_1 (u) S_2 (u) du
$$


where $S(t)$ and $\lambda(t)$ are the survival and hazard at time $t$. For the independent Cox with baseline exponential hazard (Cox-Exp) model, this becomes

$$
\pi(\tau_0, \tau_1) = \lbrack \exp \lbrace - \tau_0 (\lambda_{01} \exp (\beta'x_1)+\lambda_{02} \exp (\gamma'x_2)) \rbrace \rbrack^{-1} \int_{\tau_0}^{\tau_1} \lambda_{01} \exp (\beta'x_1)  \exp \lbrace - u (\lambda_{01} \exp (\beta'x_1)+\lambda_{02} \exp (\gamma'x_2))\rbrace du.
$$

As they are arbitrary constants, we let $\exp (\beta'x_1) = \exp (\gamma'x_2) = 1$.  Completing the integral gives a closed form expression for the absolute risk, 

$$
\pi(\tau_0, \tau_1) = \frac{\lambda_1}{\lambda_1+\lambda_2} \lbrack 1 - \exp \lbrace - (\tau_1-\tau_0)(\lambda_{01}+\lambda_{02}) \rbrace \rbrack.
$$

We define the `R` function `cox.exp` to calculate the absolute risk of disease, under the Cox-Exp model, for a specified projection.




```r
cox.exp <- function(tau0, tau1, lambda1, lambda2) {
    lambda1/(lambda1 + lambda2) * (1 - exp(-(lambda1 + lambda2) * (tau1 - tau0)))
}
```




We can perform a brief simulation to validate the expression. The function `rcox.exp` is a random bivariate number generator for the disease and competing risk event times. 




```r
rcox.exp <- function(n, lambda1, lambda2) {
    T1 <- rexp(n, rate = lambda1)
    T2 <- rexp(n, rate = lambda2)
    cbind(T1, T2)
}
```




We generate 10,000 samples from the Cox-Exp model, with the expectation for the disease event set to $E \lbrack T_1 \rbrack= 1$ and the competing event to $E \lbrack T_2 \rbrack= 2$, and obtain a Monte-Carlo estimate of the primary event occurring between the time $[1, 3)$ before the competing event occurs and given event-free survival up to time $1$.



```r
data <- rcox.exp(n = 10000, lambda1 = 1, lambda2 = 0.5)
data <- data[data[, 1] >= 1 & data[, 2] >= 1, ]
head(data)
```

```
##         T1    T2
## [1,] 1.393 2.387
## [2,] 3.145 2.122
## [3,] 2.450 1.810
## [4,] 2.108 4.411
## [5,] 1.686 1.070
## [6,] 1.269 4.363
```

```r
mean(data[, 1] < data[, 2] & data[, 1] >= 1 & data[, 1] < 3)
```

```
## [1] 0.6391
```




Now, we check this against the analytic result.



```r
cox.exp(1, 3, lambda1 = 1, lambda2 = 0.5)
```

```
## [1] 0.6335
```




Next, we generalize the random number generator and absolute risk calculation to allow for multiple competing events. The random number generator is then used to simulate the outcome of a cohort study with two competing events and three unassociated covariates, which will be used for testing the programming of the `coxph.risk` functions.



```r
cox.exp <- function(tau0, tau1, lambdas) {
    lambdas[1]/sum(lambdas) * (1 - exp(-(sum(lambdas)) * (tau1 - tau0)))
}
```




As a quick check for an error, we compute the same absolute risk in the earlier example.



```r
cox.exp(1, 3, c(1, 0.5, 0))
```

```
## [1] 0.6335
```




A similar modification is made to the random number generator.



```r
rcox.exp <- function(n, lambdas) {
    
    Ts <- sapply(lambdas, function(lambda, n) {
        rexp(n, rate = lambda)
    }, n = n)
    
    Ts
}
```




From this, we can construct a test data frame with 10,000 observations, one disease event and two competing risks. We will assume no loss to follow-up. We use `pmin` to identify the earliest time out of each row's disease and competing events. The earliest time is the only observed time for each observation (that's what makes them competing events!). It will be the observed time for one of the event types and the censoring time for all others.



```r
lambdas <- c(2, 0.3, 0.5)
times <- rcox.exp(n = 3000, lambdas)

min.time <- pmin(times[, 1], times[, 2], times[, 3])

test <- data.frame(time = min.time, disease = times[, 1] == min.time, 
    competing1 = times[, 2] == min.time, competing2 = times[, 3] == min.time, 
    x1 = rnorm(nrow(times)), x2 = runif(nrow(times)), x3 = rbinom(n = nrow(times), 
        size = 1, p = 0.3))

head(test)
```

```
##      time disease competing1 competing2      x1      x2 x3
## 1 0.04259    TRUE      FALSE      FALSE  0.5341 0.23114  0
## 2 0.29454    TRUE      FALSE      FALSE  0.3914 0.31599  0
## 3 0.40867   FALSE      FALSE       TRUE  0.4941 0.04549  0
## 4 0.37190    TRUE      FALSE      FALSE -0.7755 0.67611  0
## 5 0.15547    TRUE      FALSE      FALSE -1.1249 0.39124  0
## 6 0.13071    TRUE      FALSE      FALSE  0.1208 0.02004  0
```




***

Estimating absolute risk with `coxph.risk`
========================================================

The `coxph.risk` computes the absolute risk of an event occurring between time [time0, time1) in the presence of competing events and given event-free survival up to time time0. Although the absolute risk estimator (of section?) can apply to a more general class of relative risk models, including models with non-linear effects, the `coxph.risk` implementation is for the Cox proportional hazards model only. 

To demonstrate the usage of this function, we will compute predicted risks based on the `test` data set created in the example above. The `test` data set contains outcomes for three time-to-events, which each follow a Cox-Exp model, having the baseline hazard rates of `2`, `0.3`, and `0.5` and three unassociated (noise) covariates `x1`, `x2` and `x3`. 

After running R CMD INSTALL on the source package we load `coxph.risk`.


```r
library(coxph.risk)
```

```
## Loading required package: survival
```

```
## Loading required package: splines
```




We can verify the simulation procedure with the parametric regression function `survreg`. The exponential distribution is equivalent to a Weibull with scale parameter set to 1. For the parameterization of the `dexp` and `dweibull` distributions in R, the Weibull shape parameter is 1/rate of the exponential rate.



```r
fit1 <- survreg(Surv(time, disease) ~ x1 + x2 + x3, dist = "exponential", 
    data = test)

fit2 <- survreg(Surv(time, competing1) ~ x1 + x2 + x3, dist = "exponential", 
    data = test)

fit3 <- survreg(Surv(time, competing2) ~ x1 + x2 + x3, dist = "exponential", 
    data = test)

1/exp(coef(fit1)[1])
```

```
## (Intercept) 
##       2.078 
```

```r
1/exp(coef(fit2)[1])
```

```
## (Intercept) 
##      0.2627 
```

```r
1/exp(coef(fit3)[1])
```

```
## (Intercept) 
##      0.5088 
```




The first step in the absolute risk calculation is fitting the relative risk models. Here we suppose that the following Cox proportional hazards models were the selected relative risk models for each event type.



```r
cox1 <- coxph(Surv(time, disease) ~ x2, data = test)
cox2 <- coxph(Surv(time, competing1) ~ x1 + x2 + x3, data = test)
cox3 <- coxph(Surv(time, competing2) ~ x1, data = test)
```




Although the covariates are not associated with any of the event times, we have fit different hazard models for each to emphasize that the models need not share risk factors or in any way be dependent on each other.

To make the absolute risk calculation, in addition to the Cox models, we must specify the projection interval and the risk type for the individual for whom the projection is being made. The risk type is determined by the covariate pattern, which is supplied by a data frame containing the variables in each of the primary event and competing event models. Currently, the function supports the risk calculation for only a single risk type.

In the following, we compute several risk predictions for the risk type corresponding to the first row of the `test` data set. The primary event is indicated by the order of the Cox models, with the first position designating the primary event.

For this example, we show the syntax for estimating the absolute risk of disease between time [0, 1), for an individual with no exposures (the reference group).




```r
temp <- test[1, ]
temp[1, c("x1", "x2", "x3")] <- rep(0, 3)
coxph.risk(0, 1, projection.data = temp, cox1, cox2, cox3)
```

```
## [1] 0.6823
```




This can be compared to the analytic risk calculation for the exponential model with baseline hazard rates used to simulate the `test` data set.



```r
cox.exp(0, 1, lambdas)
```

```
## [1] 0.6708
```




By changing the period of risk we see how surviving to the beginning of the projection interval impacts risk.



```r
coxph.risk(0.2, 1, projection.data = temp, cox1, cox2, cox3)
```

```
## [1] 0.649
```

```r
coxph.risk(0.4, 1, projection.data = temp, cox1, cox2, cox3)
```

```
## [1] 0.6081
```

```r
coxph.risk(0.6, 1, projection.data = temp, cox1, cox2, cox3)
```

```
## [1] 0.5118
```






```r
cox.exp(0.2, 1, lambdas)
```

```
## [1] 0.6382
```

```r
cox.exp(0.4, 1, lambdas)
```

```
## [1] 0.5812
```

```r
cox.exp(0.6, 1, lambdas)
```

```
## [1] 0.4812
```



