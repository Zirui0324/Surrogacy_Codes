Split2 <- function(K, input) {
  I1 <- rep(0, input)        # initiate the fold labels storage for each data point
  half <- input / 2
  
  # random permutation of [1~half]
  firstHalfIDs <- sample(1:half, half, replace = FALSE)
  
  # number of points per fold in the first half
  m <- half / K
  
  # assign folds to the first half
  for (i in 1:K) {
    # The slice for fold i
    theseIDs <- firstHalfIDs[((i - 1) * m + 1) : (i * m)]
    I1[theseIDs] <- i
  }
  
  # random permutation of [half+1~input]
  secondHalfIDs <- sample((half + 1) : input, half, replace = FALSE)
  
  # assign folds to the second half
  for (i in 1:K) {
    theseIDs <- secondHalfIDs[((i - 1) * m + 1) : (i * m)]
    I1[theseIDs] <- i
  }
  
  return(I1)
}
