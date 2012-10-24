function PathLogPrior = ComputePathLogPrior(Z)

PathLogPrior = sum(log(normpdf(Z,0,1)));