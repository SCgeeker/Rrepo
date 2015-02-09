# Retrieve from http://www.agr.kuleuven.ac.be/vakken/statisticsbyR/ANOVAbyRr/multiplecompJIMRC.htm
# compiled by Jim Robison-Cox

all.pairs <- function(r)
  list(first = rep(1:r,rep(r,r))[lower.tri(diag(r))],
       second = rep(1:r, r)[lower.tri(diag(r))])

tukeyCI <- function(fitted, nis, df, MSE, conf.level=.95){
  ## fitted is a sequence of means
  ## nis is a corresponding sequence of sample sizes for each mean
  ## df is the residual df from the ANOVA table
  ## MSE = mean squared error from the ANOVA table
  ## conf.level is the family-wise confidence level, defaults to .95
  r <- length(fitted)
  pairs <- all.pairs(r)
  diffs <- fitted[pairs$first] - fitted[pairs$second]
  df <- sum(nis) - r
  T <- qtukey(conf.level, r, df)/sqrt(2)
  hwidths <-  T*sqrt(MSE*(1/nis[pairs$first] + 1/nis[pairs$second]))
  val <- cbind(diffs - hwidths, diffs, diffs + hwidths)
  dimnames(val) <- list(paste("mu",pairs$first," - mu", pairs$second,
  sep=""), c("Lower", "Diff","Upper"))
  val
}

scheffeCI <- function(fitted, nis, df, MSE, conf.level=.95){
  ## fitted is a sequence of means
  ## nis is a corresponding sequence of sample sizes for each mean
  ## df is the residual df from the ANOVA table
  ## MSE = mean squared error from the ANOVA table
  ## conf.level is the family-wise confidence level, defaults to .95
  r <- length(fitted)
  pairs <- all.pairs(r)
  diffs <- fitted[pairs$first] - fitted[pairs$second]
  T <- sqrt((r-1)*qf(conf.level,r-1,df))
  hwidths <-  T*sqrt(MSE*(1/nis[pairs$first] + 1/nis[pairs$second]))
  val <- cbind(diffs - hwidths, diffs, diffs + hwidths)
  dimnames(val) <- list(paste("mu",pairs$first," - mu", pairs$second,
  sep=""), c("Lower", "Diff","Upper"))
  val
}

bonferroniCI <- function(fitted, nis, df, MSE, conf.level=.95){
  ## fitted is a sequence of means
  ## nis is a corresponding sequence of sample sizes for each mean
  ## df is the residual df from the ANOVA table
  ## MSE = mean squared error from the ANOVA table
  ## conf.level is the family-wise confidence level, defaults to .95
  r <- length(fitted)
  pairs <- all.pairs(r)
  diffs <- fitted[pairs$first] - fitted[pairs$second]
  T <- qt(1-(1-conf.level)/(2*r*(r-1)),df)
  hwidths <-  T*sqrt(MSE*(1/nis[pairs$first] + 1/nis[pairs$second]))
  val <- cbind(diffs - hwidths, diffs, diffs + hwidths)
  dimnames(val) <- list(paste("mu",pairs$first," - mu", pairs$second,
  sep=""), c("Lower", "Diff","Upper"))
  val
}


####### DEMO
# We prepare the agrument list for the functions: 
# 
#     *
# 
#       averages,
# 
#     *
# 
#       the number of replications, 
# 
#     *
# 
#       dfMSE and 
# 
#     *
# 
#       MSE 
# 
# for Rust as function of each factor-level Brand.
# 
# > Rust.means <- tapply(Rust,Brand,mean);Rust.means
#      1      2      3      4 
#  43.14  89.44  67.95  40.47 
# 
# > Rust.len <- tapply(Rust,Brand,length);Rust.len
#        1   2   3   4 
#       10  10  10  10 
# 
# > dfMSE=fm1$df.residual;dfMSE
# [1] 36
# 
# > MSE=sum(fm1$residuals^2)/dfMSE;MSE
# [1] 6.139833