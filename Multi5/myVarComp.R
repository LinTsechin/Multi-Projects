myVarComp <- function(CovList, SampleNum) {
  if (is.list(CovList)) {
    if (length(CovList) < 2)
      stop("CovList must be a list with at least 2 elements")
    ps <- as.vector(sapply(CovList, dim))
    if (sum(ps[1] == ps) != length(ps))
      stop("All covariance matrices must have the same dimension")
    p <- ps[1]
    k <- length(CovList)
    if (length(SampleNum) < 2) {
      stop("need at least 2 populations")
    } else if (length(SampleNum) != k) {
      stop("n must be equal length(covmat)")
    }
    DNAME <- deparse(substitute(CovList))
  }
  else
    stop("covmat must be a list")
  
  n.t <- SampleNum - 1
  n <- sum(n.t)
  A.t <- lapply(1:length(CovList), 
                function(i, mat, n) { n[i] * mat[[i]] }, mat=CovList, n=n.t)
  A <- matrix(colSums(matrix(unlist(A.t), ncol=p^2, byrow=T)), ncol=p)
  logS.t <- sapply(1:length(CovList), 
                   function(i, mat, n) {n[i] * log(det(mat[[i]]))}, mat=CovList, n=n.t)
  logS <- (n-k)*log(det(A/(n-k)))
  f <- (k-1)*p*(p+1)/2
  
  if(sum(n.t[1] == n.t) < length(n.t)) {
    d <- (2*p^2+3*p-1)/(6*(p+1)*(k-1))*(sum(1/n.t)-1/(n-k))
  } else {
    d <- (2*p^2+3*p-1)*(k+1)/(6*(p+1)*(n-k))
  }
    
  STATISTIC <- (1-d)*(logS - sum(logS.t))
  PVAL <- 1 - pchisq(STATISTIC, f)
  names(STATISTIC) <- "corrected lambda*"
  names(f) <- "df"
  RVAL <- structure(list(statistic = STATISTIC, parameter = f, 
                         p.value = PVAL, data.name = DNAME, 
                         method = "Equality of Covariances Matrices Test"),
                    class="htest")
  return(RVAL)
}