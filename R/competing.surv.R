# COMPUTE INDIVIDUALIZED SURVIVAL
competing.surv <- function(baseline, times, rr=1){
  H <- matrix(baseline$haz, nrow(baseline), ncol = length(times))
  AtRisk <- outer(baseline$time, times, function(x, y) x<=y)
exp(-colSums(H*AtRisk))^rr
}