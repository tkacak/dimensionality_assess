#Empirical Kaiser Criterion
EKC <- function(data, restricted = TRUE, use = "pairwise.complete.obs"){
  p <- ncol(data)
  n <- nrow(data)
  l <- (1 + sqrt(p/n))^2
  e <- eigen(cor(data, use = use))$values
  v <- c(0, cumsum(e)[1:p - 1])
  w <- p:1
  if(restricted == TRUE){
    ref <- pmax(rep(1,p),(((p - v)/w) * l))
  }
  else {
    ref <- (((p - v)/w) * l)
  }
  which(e <= ref)[1] - 1
}
