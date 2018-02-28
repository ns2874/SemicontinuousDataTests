#' Simulated Power Calculations for Semicontinuous data
#'
#' This function runs a simulation to determine power for various tests under 
#' inputted conditions
#' @param prob.x,prob.y,mu.x,mu.y,sd.x,sd.y,n1,n2,simu.num,B
#' @keywords semicontinuous data hypothesis test
#' @return datafrome of power values
#' @export
#' @examples
#' CalErrorPower(.5, .7, 2, 3, 1, 1, 50, 50, 100000, 25, alpha=0.05)


########################################################################################
###### Compare the test for  semicontinuous data #######################################
### Compared Tests: ordinary t-test, wilcoxon test and permutated t-test, wilcoxon test
########################################################################################
CalErrorPower <- function(prob.x, prob.y, mu.x, mu.y, sd.x, sd.y, n1, n2, simu.num, B, alpha=0.05)
{
  pval.vec <- matrix(NA, nrow=simu.num, ncol=5)
  for(i in 1:simu.num)
  {
    ### generate the data
    zeroLoc1 <- rbinom(n1, size=1, prob=prob.x)
    xnorm <- rnorm(n1, mean=mu.x, sd=sd.x)
    x <- zeroLoc1*xnorm
    #x[zeroLoc1==1] <- 0
    zeroLoc2 <- rbinom(n2, size=1, prob=prob.y)
    ynorm <- rnorm(n2, mean=mu.y, sd=sd.y)
    y <- zeroLoc2*ynorm
    #y[zeroLoc2==1] <- 0
    
    n1 <- length(x)
    n2 <- length(y)
    n <- n1 + n2
    mu <- n1*n2/2
    sigma <- n1*n2*(n+1)/12
    
    ### Test 1: two sample t-test
    pval.vec[i,1] <- t.test(x, y)$p.value
    
    ### Test 2: wilcoxon test
    pval.vec[i,2] <- wilcox.test(x, y)$p.value
    
    ### permutate the sample
    dat <- c(x, y)
    temp.t.test <- t.test(x, y)$statistic
    temp.wil.test <- (wilcox.test(x, y)$statistic-mu)/sqrt(sigma)
    ### the proposed test: combine the discrete and continuous parts
    n.x <- c(sum(x==0), sum(x!=0))
    n.y <- c(sum(y==0), sum(y!=0))
    temp.propose.part1 <- chisq.test(rbind(n.x, n.y))$p.value
    temp.propose.part2 <- t.test(x[x!=0], y[y!=0])$p.value
    temp.propose.test <- -2*log(temp.propose.part1)-2*log(temp.propose.part2)
    perm.t <- rep(NA, B)
    perm.wil <- rep(NA, B)
    perm.propose <- rep(NA, B)
    for(b in 1:B)
    {
      dat.perm <- dat[sample(1:n, n, replace=F)]
      x.perm <- dat.perm[1:n1]
      y.perm <- dat.perm[(n1+1):n]
      perm.t[b] <- t.test(x.perm, y.perm)$statistic
      perm.wil[b] <- (wilcox.test(x.perm, y.perm)$statistic-mu)/sqrt(sigma)
      n.x.perm <- c(sum(x.perm==0), sum(x.perm!=0))
      n.y.perm <- c(sum(y.perm==0), sum(y.perm!=0))
      perm.propose.part1 <- chisq.test(rbind(n.x.perm, n.y.perm))$p.value
      perm.propose.part2 <- t.test(x.perm[x.perm!=0], y.perm[y.perm!=0])$p.value
      perm.propose[b] <- -2*log(perm.propose.part1)-2*log(perm.propose.part2)
    }
    ### Test 3: permutation t-test
    pval.vec[i,3] <- mean(perm.t < -abs(temp.t.test) | perm.t > abs(temp.t.test))
    ### Test 4: permutation wilcoxon test
    pval.vec[i,4] <-  mean(perm.wil < -abs(temp.wil.test) | perm.wil > abs(temp.wil.test))
    ### Test 5: the proposed test
    pval.vec[i,5] <- mean(perm.propose > temp.propose.test)
  }
  ### calculate the average p-value
  Res.pval <- apply(pval.vec<0.05, MARGIN=2, mean)
  return(Res.pval)
}



